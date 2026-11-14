"""What the family and cohort pages' tab badges count.

Mostly content a curator has not looked at yet. On the family page each count
covers the whole family, not the member selection, like the header QC badge,
and is blocking but pure: the page runs it off the event loop and draws the
number on the tab (``components/tab_badge.py``).

- **Sample QC**: unacknowledged high and medium flags of the family's own QC
  files. ``UPD_CANDIDATE`` flags are left to the UPD tab when there is one.
- **SNVs/Indels**: de novo calls, one per variant and carrier, of HIGH or
  MODERATE impact, with no validation in ``snvs.tsv``. With the ``denovo``
  step's tables, the candidates of its Wombat config that are strict (or not
  assessed by the step); without them, every config matching ``de_?novo``
  after the default exclusions (LCR, gnomAD-filtered, likely MNV). When
  Sample QC holds the de novo table back (a wrong relationship), the table is
  not read: the badge says why instead of counting.
- **STR**: genotypes of the newest catalog the pipeline flags (severity high
  or medium), one per locus and member, with no validation in ``strs.tsv``.
- **UPD**: the candidates no curator validated present or absent
  (:func:`genetics_viz.utils.upd.count_to_review`).

The cohort page counts over the individuals its table's filters leave, so its
badges are :class:`Pill` lists rebuilt in memory at each filter change; only
the STR flags and the run record are read, once, off the loop. Sample QC and
UPD count individuals and candidate chromosomes, split high and medium; STR
the open flagged genotypes, split the same way. The Individuals tab shows how
many the filters leave, and Workflow the last run's status.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Collection, Dict, Iterable, List, Mapping, Sequence, Tuple

from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)

from . import cohort_state, denovo, upd
from .sample_qc import (
    FlagRow,
    count_by_severity,
    gates_denovo,
    load_family_qc,
    severity_info,
    severity_summary,
)
from .sample_qc_acks import load_acks, split_acknowledged, unexplained_denovo_risk
from .sharding import get_family_path
from .str_eh import (
    eh_path,
    find_eh_labels,
    flagged_expr,
    get_cohort_eh_dir,
    get_family_eh_dir,
)
from .str_eh import read_eh_table_cached
from .str_validation import curated_statuses, str_validation_file
from .validation_status import UNVALIDATED, row_status
from .variant_rows import gnomad_filtered, in_lcr, likely_mnv, split_terms
from .vep import get_consequence_impact
from .wombat import (
    DENOVO_CONFIG,
    aggregate_by_variant,
    find_wombat_files,
    load_config_frame,
    variant_key,
)

logger = logging.getLogger(__name__)

#: Impacts a de novo call needs to count, most severe first.
DENOVO_IMPACTS = ("HIGH", "MODERATE")


@dataclass
class ReviewCount:
    """One tab's badge: how many, how severe the worst, and why."""

    n: int = 0
    #: ``high`` or ``medium``: the badge's colour.
    worst: str = "medium"
    tooltip: str = ""
    #: Nothing was counted because the content is held back; the badge says so.
    held: bool = False

    @property
    def shown(self) -> bool:
        return self.n > 0 or self.held


def qc_to_review(flags: Iterable[FlagRow], upd_tab: bool) -> ReviewCount:
    """The unacknowledged high and medium flags among ``flags``.

    ``flags`` are the family's active (unacknowledged) flags. With ``upd_tab``
    a ``UPD_CANDIDATE`` flag is left out: the UPD tab counts it.
    """
    flags = list(flags)
    kept = [f for f in flags if not (upd_tab and upd.is_upd_flag(f))]
    counts = count_by_severity(kept)
    n = counts.get("high", 0) + counts.get("medium", 0)
    if not n:
        return ReviewCount()
    tooltip = f"Sample QC: {severity_summary(counts)} not acknowledged"
    left_out = count_by_severity([f for f in flags if f not in kept])
    if left_out.get("high", 0) + left_out.get("medium", 0):
        tooltip += " (UPD candidates are counted on the UPD tab)"
    return ReviewCount(
        n=n, worst="high" if counts.get("high") else "medium", tooltip=tooltip
    )


def _impact(row: Dict) -> str:
    """The most severe of ``DENOVO_IMPACTS`` a call has, else ``""``.

    ``VEP_IMPACT`` when the table has it (comma-joined over transcripts after
    :func:`aggregate_by_variant`), else the impacts of its consequence terms.
    """
    raw = row.get("VEP_IMPACT")
    if raw:
        impacts = {p.strip().upper() for p in str(raw).split(",")}
    else:
        impacts = {
            get_consequence_impact(t) for t in split_terms(row.get("VEP_Consequence"))
        }
    return next((i for i in DENOVO_IMPACTS if i in impacts), "")


def denovo_to_review(data_dir: Path, family_id: str, cohort: Any = None) -> ReviewCount:
    """The family's de novo calls of HIGH or MODERATE impact never validated.

    With the denovo step's tables, only its strict candidates (and those it
    has not assessed) count; ``cohort`` lends its table, which knows
    ``RECURRENT``.
    """
    wombat_dir = get_family_path(data_dir, family_id) / "wombat"
    if not wombat_dir.is_dir():
        return ReviewCount()
    all_files = find_wombat_files(wombat_dir, family_id)
    evidence = denovo.load_family_denovo(data_dir, family_id, cohort)
    source = denovo.denovo_source(all_files, evidence) if evidence else None
    files = (
        [source]
        if source is not None
        else [wf for wf in all_files if DENOVO_CONFIG.search(wf["wombat_config"])]
    )
    if evidence is not None and source is None:
        evidence = None
    if not files:
        return ReviewCount()
    if gates_denovo(unexplained_denovo_risk(data_dir, family_id)):
        return ReviewCount(
            held=True,
            worst="high",
            tooltip="De novo calls held back: Sample QC finds a wrong relationship"
            " in this family",
        )
    validations = load_validation_map(data_dir / "validations" / "snvs.tsv", family_id)
    found: Dict[Tuple[str, str], str] = {}
    for wf in files:
        df = load_config_frame(wf["file_path"])
        if df is None:
            continue
        rows = aggregate_by_variant(df).to_dicts()
        if evidence is not None:
            for row in rows:
                row["Variant"] = variant_key(row)
            denovo.add_evidence(rows, evidence)
        for row in rows:
            if evidence is not None:
                if row.get("dnm_strict") == "no":
                    continue
            elif in_lcr(row) or gnomad_filtered(row) or likely_mnv(row):
                continue
            impact = _impact(row)
            if not impact:
                continue
            key, sample = variant_key(row), row.get("sample") or ""
            status: Dict = {}
            add_validation_status_to_row(status, validations, key, sample)
            if row_status(status) != UNVALIDATED:
                continue
            if found.get((key, sample)) != "HIGH":
                found[(key, sample)] = impact
    if not found:
        return ReviewCount()
    n_high = sum(1 for i in found.values() if i == "HIGH")
    which = "strict de novo candidate(s)" if evidence is not None else "de novo call(s)"
    after = "" if evidence is not None else ", after the LCR, gnomAD and MNV exclusions"
    return ReviewCount(
        n=len(found),
        worst="high" if n_high else "medium",
        tooltip=f"{len(found)} {which} of HIGH or MODERATE impact not yet"
        f" validated ({n_high} HIGH){after}",
    )


def str_to_review(data_dir: Path, family_id: str) -> ReviewCount:
    """The flagged genotypes of the newest catalog no curator validated."""
    directory = get_family_eh_dir(data_dir, family_id)
    labels = find_eh_labels(directory, family_id)
    if not labels:
        return ReviewCount()
    fam = read_eh_table_cached(eh_path(directory, family_id, labels[0], "genotypes"))
    if fam is None or len(fam) == 0:
        return ReviewCount()
    flagged = fam.filter(flagged_expr(fam))
    if len(flagged) == 0:
        return ReviewCount()
    statuses = curated_statuses(
        load_validation_map(str_validation_file(data_dir), family_id),
        flagged,
        add_validation_status_to_row,
        "Validation",
    )
    open_rows = [
        r
        for r in flagged.to_dicts()
        if (r.get("locus"), r.get("sample")) not in statuses
    ]
    if not open_rows:
        return ReviewCount()
    pairs = {(r.get("locus"), r.get("sample")) for r in open_rows}
    high = {
        (r.get("locus"), r.get("sample"))
        for r in open_rows
        if r.get("severity") == "high"
    }
    return ReviewCount(
        n=len(pairs),
        worst="high" if high else "medium",
        tooltip=f"{len(pairs)} flagged STR genotype(s) not yet validated"
        f" ({len(high)} in the pathogenic range)",
    )


def upd_to_review(data_dir: Path, family_id: str) -> ReviewCount:
    """The UPD candidates no curator validated present or absent."""
    n = upd.count_to_review(upd.family_candidates(data_dir, family_id))
    if not n:
        return ReviewCount()
    return ReviewCount(
        n=n, worst="medium", tooltip=f"{n} UPD candidate(s) no curator settled"
    )


def qc_family_to_review(data_dir: Path, family_id: str, upd_tab: bool) -> ReviewCount:
    """:func:`qc_to_review` of the family's own QC files, acknowledgements applied."""
    try:
        qc = load_family_qc(data_dir, family_id)
    except Exception:
        logger.warning("Failed to read sample QC of %s", family_id, exc_info=True)
        return ReviewCount()
    if qc is None:
        return ReviewCount()
    active, _ = split_acknowledged(qc.flags, load_acks(data_dir))
    return qc_to_review(active, upd_tab)


# ---------------------------------------------------------------------------
# Cohort page: counts over the individuals the table's filters leave
# ---------------------------------------------------------------------------

#: ``(family_id, sample)``.
Key = Tuple[str, str]

#: The severities a cohort badge counts, most severe first.
COUNTED_SEVERITIES = ("high", "medium")

#: The Individuals tab's badge colour (a Quasar colour).
FILTER_COLOR = "blue-8"


@dataclass(frozen=True)
class Pill:
    """One badge on a cohort tab: a count, or an icon, in a colour."""

    text: str
    color: str
    tooltip: str = ""
    icon: str = ""


def _plural(n: int, one: str, many: str) -> str:
    return f"{n:,} {one if n == 1 else many}"


def _severity_pills(
    counts: Mapping[str, int], tooltips: Mapping[str, str]
) -> List[Pill]:
    return [
        Pill(f"{counts[sev]:,}", severity_info(sev)["color"], tooltips[sev])
        for sev in COUNTED_SEVERITIES
        if counts.get(sev)
    ]


def filter_pills(shown: int, total: int, n_filters: int) -> List[Pill]:
    """The Individuals tab's badge: how many individuals the filters leave."""
    if not n_filters:
        return []
    return [
        Pill(
            f"{shown:,}",
            FILTER_COLOR,
            f"{shown:,} of {_plural(total, 'individual', 'individuals')} shown"
            f" ({_plural(n_filters, 'filter', 'filters')}). Every tab but"
            " Workflow shows these only.",
        )
    ]


def qc_worst_by_individual(flags: Iterable[FlagRow], upd_tab: bool) -> Dict[Key, str]:
    """Each individual's worst high or medium flag, by ``(FID, sample)``.

    ``flags`` are the cohort's active (unacknowledged) flags; a pair flag
    counts for both its samples. With ``upd_tab`` a ``UPD_CANDIDATE`` flag is
    left to the UPD tab, as on the family page.
    """
    worst: Dict[Key, str] = {}
    for f in flags:
        if upd_tab and upd.is_upd_flag(f):
            continue
        sev = severity_info(f.severity)["key"]
        if sev not in COUNTED_SEVERITIES:
            continue
        for key in f.keys:
            if worst.get(key) != "high":
                worst[key] = sev
    return worst


def qc_cohort_pills(
    worst: Mapping[Key, str], shown: Collection[Key], upd_tab: bool
) -> List[Pill]:
    """Sample QC: the shown individuals by worst unacknowledged flag."""
    counts = {sev: 0 for sev in COUNTED_SEVERITIES}
    for key in shown:
        sev = worst.get(key)
        if sev:
            counts[sev] += 1
    among = f", among the {len(shown):,} shown"
    after = " UPD candidates are counted on the UPD tab." if upd_tab else ""
    return _severity_pills(
        counts,
        {
            "high": f"{_plural(counts['high'], 'individual', 'individuals')} with"
            f" a high-severity Sample QC flag not acknowledged{among}.{after}",
            "medium": f"{_plural(counts['medium'], 'individual', 'individuals')}"
            f" whose worst unacknowledged Sample QC flag is medium{among}.{after}",
        },
    )


def upd_cohort_pills(
    open_confidences: Iterable[Sequence[str]], n_shown: int
) -> List[Pill]:
    """UPD: the candidate chromosomes no curator settled, by confidence.

    ``open_confidences`` holds, for each shown individual, the confidence of
    each of its chromosomes with an event of the shown confidences and no
    validation present or absent (the table rows' ``_upd_open``).
    """
    counts = {sev: 0 for sev in COUNTED_SEVERITIES}
    for confidences in open_confidences:
        for conf in confidences:
            counts["high" if conf == "high" else "medium"] += 1
    among = f", among the {n_shown:,} individuals shown"
    return _severity_pills(
        counts,
        {
            "high": f"{_plural(counts['high'], 'high-confidence UPD candidate', 'high-confidence UPD candidates')}"
            f" no curator settled{among}",
            "medium": f"{_plural(counts['medium'], 'medium-confidence UPD candidate', 'medium-confidence UPD candidates')}"
            f" no curator settled{among}",
        },
    )


def str_cohort_open_flags(cohort: Any, data_dir: Path) -> List[Tuple[Key, str, str]]:
    """The newest catalog's flagged genotypes no curator validated.

    ``((FID, sample), locus, severity)``, severity ``high`` or ``medium``, one
    per locus and individual. Reads the cohort's ``flags.tsv`` (a few KB),
    else its genotype table. Blocking: run it in a thread.
    """
    directory = get_cohort_eh_dir(cohort)
    labels = find_eh_labels(directory, cohort.name)
    if not labels:
        return []
    df = None
    for kind in ("flags", "genotypes"):
        df = read_eh_table_cached(eh_path(directory, cohort.name, labels[0], kind))
        if df is not None:
            break
    if df is None or len(df) == 0:
        return []
    flagged = df.filter(flagged_expr(df))
    if len(flagged) == 0:
        return []
    statuses = curated_statuses(
        load_validation_map(str_validation_file(data_dir)),
        flagged,
        add_validation_status_to_row,
        "Validation",
    )
    found: Dict[Tuple[Key, str], str] = {}
    for r in flagged.to_dicts():
        locus, sample = r.get("locus") or "", r.get("sample") or ""
        if (locus, sample) in statuses:
            continue
        ident = ((r.get("family_id") or "", sample), locus)
        if found.get(ident) != "high":
            found[ident] = "high" if r.get("severity") == "high" else "medium"
    return [(key, locus, sev) for (key, locus), sev in found.items()]


def str_cohort_pills(
    open_flags: Iterable[Tuple[Key, str, str]], shown: Collection[Key]
) -> List[Pill]:
    """STR: the shown individuals' flagged genotypes no curator validated."""
    keys = set(shown)
    counts = {sev: 0 for sev in COUNTED_SEVERITIES}
    for key, _locus, sev in open_flags:
        if key in keys:
            counts[sev] += 1
    among = f", among the {len(keys):,} individuals shown"
    return _severity_pills(
        counts,
        {
            "high": f"{_plural(counts['high'], 'flagged STR genotype', 'flagged STR genotypes')}"
            f" in the pathogenic range not yet validated{among}",
            "medium": f"{_plural(counts['medium'], 'flagged STR genotype', 'flagged STR genotypes')}"
            f" (premutation range) not yet validated{among}",
        },
    )


def denovo_cohort_status(
    cohort: Any, data_dir: Path, with_curation: bool
) -> Tuple[Dict[Key, str], List[Key]]:
    """Blocking: each child's de novo QC verdict, and the curation issues.

    The issues (``with_curation`` only) are the curated de novo calls the
    strict set hides or the candidates lack (``denovo.curation_check``), one
    ``(family, sample)`` each.
    """
    data = denovo.load_cohort_denovo(cohort)
    if data is None:
        return {}, []
    verdicts = {
        (str(q.get("family_id")), str(q.get("child"))): str(q.get("qc") or "")
        for q in data.qc_rows()
    }
    if not with_curation or not len(data.dnm):
        return verdicts, []
    families = {str(f) for f in data.dnm["family_id"].to_list()}
    curation = denovo.load_curation(data_dir, families)
    issues = [
        (r["FID"], r["sample"])
        for r in denovo.curation_check(data.rows(), curation, families)
    ]
    return verdicts, issues


def denovo_cohort_pills(
    verdicts: Mapping[Key, str], issues: Iterable[Key], shown: Collection[Key]
) -> List[Pill]:
    """De novo: children shown out of the QC range, then curation issues."""
    keys = set(shown)
    high = sum(1 for k, v in verdicts.items() if k in keys and v == "high")
    low = sum(
        1 for k, v in verdicts.items() if k in keys and v in ("low", "no_callable")
    )
    pills = []
    if high:
        pills.append(
            Pill(
                f"{high:,}",
                severity_info("high")["color"],
                f"{_plural(high, 'child', 'children')} with too many strict de novo"
                " candidates for the callable genome (QC high)",
            )
        )
    if low:
        pills.append(
            Pill(
                f"{low:,}",
                severity_info("medium")["color"],
                f"{_plural(low, 'child', 'children')} with too few strict de novo"
                " candidates, or no callable genome (QC low)",
            )
        )
    n_issues = sum(1 for k in issues if k in keys)
    if n_issues:
        pills.append(
            Pill(
                f"{n_issues:,}",
                "#6a1b9a",
                f"{_plural(n_issues, 'curated de novo call', 'curated de novo calls')}"
                " not strict, or no longer candidates",
                icon="visibility_off",
            )
        )
    return pills


def workflow_pills(
    state: "cohort_state.CohortState | None", caption: str = ""
) -> List[Pill]:
    """Workflow: the last run's status, as its chip's icon and colour.

    Nothing without a run record, or for a status this build does not know:
    the tab says it then.
    """
    run = state.last_run if state is not None and not state.unsupported else None
    if run is None or run.status not in cohort_state.RUN_STATUS_META:
        return []
    meta = cohort_state.record_status_meta(run)
    tooltip = f"Last run: {meta.get('label', run.status)}"
    if caption:
        tooltip += f" ({caption})"
    return [Pill("", meta.get("color", "grey"), tooltip, icon=meta.get("icon", "help"))]
