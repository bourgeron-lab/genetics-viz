"""The cohort page's Families tab: one row per family, without UI.

Everything a curator scans to pick the next family: structure, complete
trios, affected members, family ancestry, the highest FROH grade, Sample QC,
the de novo QC verdict, what waits for review, diagnostics, the conclusion.
The pedigree columns are computed at page build; the ancestry, ROH and de
novo tables are read off the event loop (:func:`load_family_facts`).
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set, Tuple

import polars as pl

from . import denovo
from .ancestry import ID_COLUMNS, NULL_VALUES, find_cohort_ancestry_files
from .pedigree_labels import normalize_pheno
from .report_summaries import FamilySummary, read_family_summaries
from .roh import GRADES, find_cohort_roh_files, read_roh_table
from .tsv import read_tsv_or_none
from .variant_rows import split_terms
from .vep import get_consequence_impact

logger = logging.getLogger(__name__)

Key = Tuple[str, str]

#: Impacts a strict de novo candidate needs to wait for review.
DENOVO_IMPACTS = ("HIGH", "MODERATE")


@dataclass
class FamilyFacts:
    """What the cohort tables say per individual, read once off the loop."""

    #: ``(FID, IID)`` -> predicted ancestry region.
    regions: Dict[Key, str] = field(default_factory=dict)
    #: ``(FID, IID)`` -> ``froh_category``.
    froh: Dict[Key, str] = field(default_factory=dict)
    #: ``(FID, child)`` -> de novo QC verdict (``ok``, ``high``, ``low``…).
    dn_verdicts: Dict[Key, str] = field(default_factory=dict)
    #: FID -> strict HIGH/MODERATE de novo candidates not validated,
    #: ``(n, n_high)``. Curators only.
    dn_open: Dict[str, Tuple[int, int]] = field(default_factory=dict)
    #: FID -> standing curator conclusion. Curators only.
    summaries: Dict[str, FamilySummary] = field(default_factory=dict)
    #: Which tables were found: ``ancestry``, ``roh``, ``denovo``.
    found: Set[str] = field(default_factory=set)


def _key(row: Mapping[str, Any], sample_col: str) -> Key:
    return (str(row.get("family_id") or ""), str(row.get(sample_col) or ""))


def _impact(consequences: Any) -> str:
    impacts = {get_consequence_impact(t) for t in split_terms(consequences)}
    return next((i for i in DENOVO_IMPACTS if i in impacts), "")


def dn_open_counts(
    rows: Iterable[Mapping[str, Any]], validated: Set[Key]
) -> Dict[str, Tuple[int, int]]:
    """Per family, its strict HIGH/MODERATE candidates nobody validated.

    ``validated`` holds the ``(chr:pos:ref:alt, child)`` pairs with a
    validation of any status. Counted from the cohort ``dnm`` table, so a
    Wombat row the step did not assess is not counted (the family tab's badge
    counts it).
    """
    found: Dict[str, Dict[Key, str]] = {}
    for r in rows:
        if str(r.get("strict") or "") != "yes":
            continue
        impact = _impact(r.get("consequences"))
        if not impact:
            continue
        variant = f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}:{r.get('alt')}"
        child = str(r.get("child") or "")
        if (variant, child) in validated:
            continue
        mine = found.setdefault(str(r.get("family_id") or ""), {})
        if mine.get((variant, child)) != "HIGH":
            mine[(variant, child)] = impact
    return {
        fid: (len(calls), sum(1 for i in calls.values() if i == "HIGH"))
        for fid, calls in found.items()
    }


def _validated_snvs(data_dir: Path) -> Set[Key]:
    from genetics_viz.components.validation_loader import load_validation_map

    path = data_dir / "validations" / "snvs.tsv"
    if not path.exists():
        return set()
    # An ignored record does not count: entries are (status, inheritance,
    # comment, ignore, ...).
    return {
        key
        for key, entries in load_validation_map(path).items()
        if any(entry[3] != "1" for entry in entries)
    }


def load_family_facts(cohort: Any, data_dir: Path, curator: bool) -> FamilyFacts:
    """Blocking: the cohort's ancestry, FROH and de novo tables, per individual.

    Each table is optional; one that fails to read is left out and logged.
    """
    facts = FamilyFacts()
    try:
        files = find_cohort_ancestry_files(cohort)
        path = files.ancestry_path if files else None
        df = read_tsv_or_none(path, null_values=NULL_VALUES, string_columns=ID_COLUMNS)
        if df is not None and {"IID", "family_id", "region"} <= set(df.columns):
            facts.found.add("ancestry")
            for r in df.select("family_id", "IID", "region").to_dicts():
                if r["region"]:
                    facts.regions[_key(r, "IID")] = str(r["region"])
    except Exception:
        logger.warning("Failed to read %s's ancestry", cohort.name, exc_info=True)
    try:
        files = find_cohort_roh_files(cohort)
        df = read_roh_table(files.path("froh")) if files else None
        if df is not None and {"IID", "family_id", "froh_category"} <= set(df.columns):
            facts.found.add("roh")
            for r in df.select("family_id", "IID", "froh_category").to_dicts():
                if r["froh_category"]:
                    facts.froh[_key(r, "IID")] = str(r["froh_category"])
    except Exception:
        logger.warning("Failed to read %s's FROH", cohort.name, exc_info=True)
    try:
        data = denovo.load_cohort_denovo(cohort)
        if data is not None:
            facts.found.add("denovo")
            for q in data.qc_rows():
                if q.get("qc"):
                    facts.dn_verdicts[_key(q, "child")] = str(q["qc"])
            if curator and len(data.dnm) and "strict" in data.dnm.columns:
                facts.dn_open = dn_open_counts(
                    data.dnm.select(
                        pl.col(c)
                        for c in ("family_id", "child", "chrom", "pos", "ref", "alt")
                        + ("strict", "consequences")
                        if c in data.dnm.columns
                    ).to_dicts(),
                    _validated_snvs(data_dir),
                )
    except Exception:
        logger.warning("Failed to read %s's de novo tables", cohort.name, exc_info=True)
    if curator:
        try:
            facts.summaries = read_family_summaries(data_dir)
        except Exception:
            logger.warning("Failed to read the family conclusions", exc_info=True)
    return facts


def grade_rank(category: Optional[str]) -> int:
    """A ``froh_category``'s strength: the strongest grade highest, none 0."""
    order = [g.category for g in GRADES]
    return len(order) - order.index(category) if category in order else 0


def is_affected(phenotype: Optional[str]) -> bool:
    return normalize_pheno(phenotype) == "2"


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------

#: Structure kinds, in the filter's order.
STRUCTURES = ("Singleton", "Duo", "Trio", "Quartet", "Extended", "No parents")
#: The filter options of a column with no value for a family.
NOT_ASSESSED = "Not assessed"
NO_DIAGNOSTIC = "None"
#: Diagnostic verdicts, strongest first, and their colours (the individuals').
DIAG_ORDER = ("pathogenic", "uncertain", "benign")
DIAG_COLORS = {"pathogenic": "#dc2626", "uncertain": "#f59e0b", "benign": "#16a34a"}
#: De novo QC verdicts, worst first.
DN_ORDER = ("high", "low", "no_callable", "ok")


def family_structure(family: Any) -> Tuple[str, str, int]:
    """``(kind, label, complete trios)`` from the sequenced members.

    A complete trio is a child with both parents sequenced. Dummies (never
    sequenced) do not count, so twins with dummy parents are "No parents".
    """
    members = list(family.samples)
    ids = {s.sample_id for s in members}
    n = len(members)
    trio_children = [s for s in members if s.father_id in ids and s.mother_id in ids]
    with_parent = [s for s in members if s.father_id in ids or s.mother_id in ids]
    trios = len(trio_children)
    if n == 1:
        kind = "Singleton"
    elif not with_parent:
        kind = "No parents"
    elif n == 2:
        kind = "Duo"
    elif n == 3 and trios == 1:
        kind = "Trio"
    elif (
        n == 4
        and trios == 2
        and len({(s.father_id, s.mother_id) for s in trio_children}) == 1
    ):
        kind = "Quartet"
    else:
        kind = "Extended"
    label = f"{kind} ({n})" if kind in ("Extended", "No parents") else kind
    return kind, label, trios


def family_ancestry(family: Any, regions: Mapping[Key, str]) -> str:
    """The founders' regions, most common first (``EUR / AFR``).

    Founders (no sequenced parent) say where the family comes from; with none
    assigned, every member's region is used.
    """
    fid = family.family_id
    founders = [s for s in family.samples if s.is_founder]
    found = [regions.get((fid, s.sample_id)) for s in founders]
    if not any(found):
        found = [regions.get((fid, s.sample_id)) for s in family.samples]
    counts: Dict[str, int] = {}
    for region in found:
        if region:
            counts[region] = counts.get(region, 0) + 1
    return " / ".join(sorted(counts, key=lambda r: (-counts[r], r)))


def _plural(n: int, one: str, many: str) -> str:
    return f"{n:,} {one if n == 1 else many}"


def family_row(
    family: Any,
    facts: FamilyFacts,
    *,
    qc: Optional[Dict[str, Any]] = None,
    upd_open: Iterable[str] = (),
    str_open: Tuple[int, int] = (0, 0),
    diagnostics: Tuple[int, str] = (0, ""),
    curator: bool = False,
    show_relationship: bool = False,
) -> Dict[str, Any]:
    """The Families table's row of one family.

    ``qc`` is the family's Sample QC cell (``QC``, ``QC_badges``, ``_qc_cats``,
    ``_qc_rank``), built by the page as the individuals' are. ``upd_open``
    lists the confidence of each unsettled UPD candidate, ``str_open`` is
    ``(n, n_high)`` of open STR flags and ``diagnostics`` ``(count, strongest
    verdict)``. Review and conclusion fields are filled for curators only.
    """
    from .roh import grade_badge, grade_for
    from .sample_qc import severity_info

    fid = family.family_id
    kind, label, trios = family_structure(family)
    members = list(family.samples)
    affected = sum(1 for s in members if is_affected(s.phenotype))
    row: Dict[str, Any] = {
        "FID": fid,
        "Structure": label,
        "_structure": kind,
        "Members": len(members),
        "Trios": trios,
        "Affected": f"{affected} / {len(members)}",
        "_affected": affected,
        "Ancestry": family_ancestry(family, facts.regions),
        **(qc or {}),
    }

    # The strongest FROH grade of any member.
    categories = [facts.froh.get((fid, s.sample_id)) for s in members]
    best = max(categories, key=grade_rank, default=None)
    if "roh" in facts.found and best is not None:
        row["FROH"] = ""
        row["FROH_badges"] = [grade_badge(best, show_relationship)]
        row["_froh_rank"] = grade_rank(best)
        row["_froh_grade"] = grade_for(best).grade
    else:
        row.update(FROH="", FROH_badges=[], _froh_rank=-1, _froh_grade=NOT_ASSESSED)

    # The worst de novo QC verdict among the children.
    verdicts = {
        s.sample_id: facts.dn_verdicts[(fid, s.sample_id)]
        for s in members
        if (fid, s.sample_id) in facts.dn_verdicts
    }
    if verdicts:
        worst = min(
            verdicts.values(),
            key=lambda v: DN_ORDER.index(v) if v in DN_ORDER else len(DN_ORDER),
        )
        info = denovo.qc_info(worst)
        row["DenovoQC"] = ""
        row["DenovoQC_badges"] = [
            {
                "label": info["label"],
                "color": info["color"],
                "tooltip": "\n".join(
                    f"{child}: {denovo.qc_info(v)['label']}"
                    for child, v in sorted(verdicts.items())
                ),
            }
        ]
        row["_dn_rank"] = len(DN_ORDER) - (
            DN_ORDER.index(worst) if worst in DN_ORDER else len(DN_ORDER)
        )
        row["_dn_cats"] = sorted(
            {denovo.qc_info(v)["label"] for v in verdicts.values()}
        )
    else:
        row.update(
            DenovoQC="", DenovoQC_badges=[], _dn_rank=-1, _dn_cats=[NOT_ASSESSED]
        )

    # Standing diagnostics: how many, and the strongest verdict.
    n_diag, strongest = diagnostics
    if n_diag:
        row["Diagnostics"] = ""
        row["Diagnostics_badges"] = [
            {
                "label": f"{strongest.capitalize()} · {n_diag}"
                if n_diag > 1
                else strongest.capitalize(),
                "color": DIAG_COLORS.get(strongest, "#9e9e9e"),
                "tooltip": f"{_plural(n_diag, 'diagnostic', 'diagnostics')};"
                f" the strongest verdict is {strongest}",
            }
        ]
        row["_diag_rank"] = len(DIAG_ORDER) - (
            DIAG_ORDER.index(strongest) if strongest in DIAG_ORDER else len(DIAG_ORDER)
        )
        row["_diag_cats"] = [strongest]
    else:
        row.update(Diagnostics="", Diagnostics_badges=[], _diag_rank=0)
        row["_diag_cats"] = [NO_DIAGNOSTIC]

    if not curator:
        return row

    # What waits for a curator: the family tabs' badges, from cohort tables.
    review: List[Dict[str, str]] = []
    total = 0
    dn_n, dn_high = facts.dn_open.get(fid, (0, 0))
    upd_list = list(upd_open)
    upd_high = sum(1 for c in upd_list if c == "high")
    str_n, str_high = str_open
    for name, n, n_high, what in (
        ("De novo", dn_n, dn_high, "strict de novo candidates of HIGH or MODERATE"
         " impact not validated"),
        ("UPD", len(upd_list), upd_high, "UPD candidates no curator settled"),
        ("STR", str_n, str_high, "flagged STR genotypes not validated"),
    ):  # fmt: skip
        if not n:
            continue
        total += n
        review.append(
            {
                "label": f"{name} {n}",
                "color": severity_info("high" if n_high else "medium")["color"],
                "tooltip": f"{n} {what}",
            }
        )
    row["Review"] = ""
    row["Review_badges"] = review
    row["_review_n"] = total

    summary = facts.summaries.get(fid)
    row["Conclusion"] = summary.date if summary else ""
    row["Conclusion_badges"] = (
        [
            {
                "label": summary.date,
                "color": "#546e7a",
                "tooltip": f"{summary.user}, {summary.date}:\n{summary.text[:300]}",
            }
        ]
        if summary
        else []
    )
    return row


# ---------------------------------------------------------------------------
# Columns and filters
# ---------------------------------------------------------------------------


def _badges(column_id: str, header: str, sort_field: str, **extra: Any) -> Dict:
    return {
        "id": column_id,
        "header": header,
        "cellType": "badge_list",
        "sortable": True,
        "sortField": sort_field,
        "sorting": "numerical",
        **extra,
    }


def _multi(options: Iterable[str]) -> Dict[str, Any]:
    return {"type": "multiselect", "options": list(options), "placeholder": "All"}


def family_columns(
    cohort_name: str,
    *,
    qc_options: Optional[List[str]] = None,
    diag_options: Iterable[str] = (),
    with_ancestry: bool = False,
    with_roh: bool = False,
    with_denovo: bool = False,
    curator: bool = False,
) -> List[Dict[str, Any]]:
    """The Families table's columns; an analysis column only when it ran.

    ``qc_options`` are the individuals table's Sample QC filter options (none
    when the cohort has no Sample QC).
    """
    from .roh import GRADES

    columns: List[Dict[str, Any]] = [
        {
            "id": "FID",
            "header": "Family ID",
            "cellType": "link",
            "href": f"/cohort/{cohort_name}/family/{{FID}}",
            "sortable": True,
            "minWidth": 90,
            "filter": {"type": "text", "placeholder": "Filter...", "size": 8},
        },
        {
            "id": "Structure",
            "header": "Structure",
            "sortable": True,
            "filter": _multi(STRUCTURES),
        },
        {
            "id": "Trios",
            "header": "Trios",
            "sortable": True,
            "sorting": "numerical",
            "filter": {"type": "checkbox", "label": "Has a trio"},
        },
        {
            "id": "Affected",
            "header": "Affected",
            "sortable": True,
            "sortField": "_affected",
            "sorting": "numerical",
            "filter": {"type": "checkbox", "label": "2+ affected"},
        },
    ]
    if with_ancestry:
        columns.append(
            {
                "id": "Ancestry",
                "header": "Ancestry",
                "sortable": True,
                "filter": {"type": "text", "placeholder": "Filter...", "size": 6},
            }
        )
    if with_roh:
        grades = list(dict.fromkeys(g.grade for g in GRADES))
        columns.append(
            _badges(
                "FROH", "FROH", "_froh_rank", filter=_multi(grades + [NOT_ASSESSED])
            )
        )
    if qc_options:
        columns.append(
            _badges("QC", "Sample QC", "_qc_rank", filter=_multi(qc_options))
        )
    if with_denovo:
        labels = list(dict.fromkeys(denovo.qc_info(v)["label"] for v in DN_ORDER))
        columns.append(
            _badges(
                "DenovoQC",
                "De novo QC",
                "_dn_rank",
                filter=_multi(labels + [NOT_ASSESSED]),
            )
        )
    if curator:
        columns.append(
            _badges(
                "Review",
                "To review",
                "_review_n",
                filter={"type": "checkbox", "label": "Waiting"},
            )
        )
    columns.append(
        _badges(
            "Diagnostics",
            "Diagnostics",
            "_diag_rank",
            filter=_multi(list(diag_options) + [NO_DIAGNOSTIC]),
        )
    )
    if curator:
        columns.append(
            {
                "id": "Conclusion",
                "header": "Conclusion",
                "cellType": "badge_list",
                "sortable": True,
                "filter": {"type": "checkbox", "label": "None yet"},
            }
        )
    return columns


def filter_families(
    rows: List[Dict[str, Any]], filters: Mapping[str, Any], qc_any: str = ""
) -> List[Dict[str, Any]]:
    """The families the table's column filters leave.

    ``qc_any`` is the Sample QC option meaning "any unacknowledged flag".
    """
    out = rows
    for col in ("FID", "Ancestry"):
        text = (filters.get(col) or "").strip().lower()
        if text:
            out = [r for r in out if text in str(r.get(col) or "").lower()]
    for col, field_ in (
        ("Structure", "_structure"),
        ("FROH", "_froh_grade"),
    ):
        wanted = set(filters.get(col) or [])
        if wanted:
            out = [r for r in out if r.get(field_) in wanted]
    for col, field_ in (("DenovoQC", "_dn_cats"), ("Diagnostics", "_diag_cats")):
        wanted = set(filters.get(col) or [])
        if wanted:
            out = [r for r in out if wanted & set(r.get(field_) or [])]
    wanted_qc = set(filters.get("QC") or [])
    if wanted_qc:
        out = [
            r
            for r in out
            if wanted_qc & set(r.get("_qc_cats") or [])
            or (qc_any in wanted_qc and r.get("_qc_rank", -1) > 0)
        ]
    if filters.get("Trios"):
        out = [r for r in out if r.get("Trios", 0) > 0]
    if filters.get("Affected"):
        out = [r for r in out if r.get("_affected", 0) >= 2]
    if filters.get("Review"):
        out = [r for r in out if r.get("_review_n", 0) > 0]
    if filters.get("Conclusion"):
        out = [r for r in out if not r.get("Conclusion")]
    return out


#: The genome QC flag that marks a sample's DNA as from a cell line.
CELL_LINE = "CELL_LINE"


def compact_qc_cell(
    cell: Dict[str, Any],
    flags: List[Any],
    acked: List[Any],
    family_id: str,
    acked_label: str,
) -> Dict[str, Any]:
    """A family's Sample QC cell as one badge, plus a cell-line badge.

    ``cell`` is the cell built as an individual's (one badge per category);
    its categories stay in ``_qc_cats`` for the filter, and in the tooltip.
    The badge reads the worst severity and how many members it flags. A
    ``CELL_LINE`` flag, acknowledged or not, keeps a badge of its own: it is a
    caveat on the family's results, not only an issue to review.
    """
    from .sample_qc import severity_info

    badges = cell.get("QC_badges") or []
    shown: List[Dict[str, str]] = []
    if flags:
        worst = severity_info(max(flags, key=lambda f: f.rank).severity)
        members = {s for f in flags for fid, s in f.keys if fid == family_id and s}
        active = [b for b in badges if b.get("label") != acked_label]
        shown.append(
            {
                "label": f"{worst['label']} · {_plural(len(members), 'member', 'members')}",
                "color": worst["color"],
                "tooltip": "\n".join(
                    f"{b['label']}: {b.get('tooltip', '')}".rstrip(": ") for b in active
                ),
            }
        )
    elif acked:
        shown.extend(b for b in badges if b.get("label") == acked_label)
    if any(f.code == CELL_LINE for f in list(flags) + list(acked)):
        shown.append(
            {
                "label": "Cell line",
                "color": "#6a1b9a",
                "tooltip": "DNA from a cell line: de novo, UPD and CNV results"
                " need confirmation on primary DNA",
            }
        )
    return {**cell, "QC_badges": shown}


# ---------------------------------------------------------------------------
# Previous / next on the family page
# ---------------------------------------------------------------------------

#: ``app.storage.user`` key: per cohort, the Families table's family order
#: when it is filtered or sorted. The family page's arrows walk it.
FAMILY_ORDER_KEY = "family_order"


def save_family_order(
    storage: Any, cohort_name: str, order: Optional[List[str]]
) -> None:
    """Remember ``order`` for ``cohort_name``; ``None`` forgets it.

    A new dict is assigned, so a persisted storage sees the change.
    """
    orders = dict(storage.get(FAMILY_ORDER_KEY) or {})
    if order is None:
        if cohort_name not in orders:
            return
        orders.pop(cohort_name)
    elif orders.get(cohort_name) == order:
        return
    else:
        orders[cohort_name] = list(order)
    storage[FAMILY_ORDER_KEY] = orders


def family_order(
    storage: Any, cohort_name: str, family_id: str, all_families: List[str]
) -> Tuple[List[str], bool]:
    """The families the arrows walk, and whether they are the Families tab's.

    The remembered order when it holds ``family_id``; else every family by
    ID (a family opened from Search, a gene page or a shared link).
    """
    saved = (storage.get(FAMILY_ORDER_KEY) or {}).get(cohort_name) or []
    if family_id in saved:
        return list(saved), True
    return sorted(all_families), False
