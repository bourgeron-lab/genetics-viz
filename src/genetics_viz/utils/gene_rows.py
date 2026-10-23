"""Gene and variant pages without UI: index rows to table rows, curation joined.

The rows come from :mod:`genetics_viz.utils.gene_index` (every cohort of every
data directory). Validations and diagnostics are read live from each data
directory's files and joined by the same keys the family tabs save under, so a
curator's verdict on a family page shows here at once.

A row of another data directory than the user's own links through
``/goto/<n>``, which switches the user's data directory before opening the
family: a family page reads the user's current one.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import polars as pl

from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.utils.data import (
    get_all_data_stores,
    get_data_dir_options,
    get_static_prefix_map,
)
from genetics_viz.utils.diagnostic_files import diagnostic_path, diagnostic_paths
from genetics_viz.utils.gene_scoring import build_gene_badges, split_genes
from genetics_viz.utils.str_eh import severity_info, status_info
from genetics_viz.utils.validation_files import validation_path
from genetics_viz.utils.validation_status import UNVALIDATED
from genetics_viz.utils.variant_rows import add_vep_badges, consequence_badges

#: Curation kinds joined to the index kinds (UPDs name no gene).
_CURATED = {"snv": "SNV", "sv": "SV", "str": "STR"}

#: The family page tab each kind opens on.
_FAMILY_TAB = {"snv": "wombat", "sv": "svs", "str": "str", "roh": "roh"}


# ---------------------------------------------------------------------------
# Data directories
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DirInfo:
    """A registered data directory, as the rows name and link it."""

    path: str
    #: N of its ``/data-N`` prefix, used by ``/goto/<n>``: a link carries no path.
    index: str
    label: str


def data_dir_infos() -> Dict[str, DirInfo]:
    """Every registered data directory, by path."""
    labels = {o["value"]: o["label"] for o in get_data_dir_options()}
    return {
        path: DirInfo(
            path, prefix.rsplit("-", 1)[-1], labels.get(path, Path(path).name)
        )
        for path, prefix in get_static_prefix_map().items()
    }


def goto_fields(row: Dict[str, Any], info: Optional[DirInfo], next_path: str) -> None:
    """``_dir`` and ``_next`` for a ``/goto/{_dir}?next={_next}`` link cell,
    and the row's ``_key``: the same family link can recur in a table."""
    row["_dir"] = info.index if info else "0"
    row["_next"] = next_path
    row["_key"] = (
        "|".join(
            str(row.get(k) or "")
            for k in (
                "data_dir",
                "Type",
                "Cohort",
                "FID",
                "Sample",
                "Variant",
                "Segment",
            )
        )
        + f"|{row.get('Timestamp') or ''}"
    )


def gene_keys(symbol: str, snapshot: Any) -> Tuple[str, Set[str]]:
    """The gene's current symbol, and every upper-cased name it goes by.

    VEP symbols lag HGNC, so a table may hold ``KIAA1109`` for ``BLTP1``:
    previous and alias symbols are looked up too, when they name one gene
    (``GenetrekSnapshot.matching_identifiers``). Without a snapshot, the
    symbol as typed.
    """
    typed = symbol.strip()
    current = snapshot.resolve(typed) or typed.upper()
    keys = {typed.upper(), current.upper()}
    keys |= set(snapshot.matching_identifiers({current}))
    return current, keys


def family_path(cohort: str, family_id: str, kind: str, config: str = "") -> str:
    """The family page, on the tab of ``kind`` (and the Wombat ``config``)."""
    path = f"/cohort/{cohort}/family/{family_id}?tab={_FAMILY_TAB[kind]}"
    if kind == "snv" and config:
        path += f"&wcfg={config}"
    return path


# ---------------------------------------------------------------------------
# Curation
# ---------------------------------------------------------------------------


@dataclass
class Curation:
    """Validation and diagnostic lookups of every data directory, per kind."""

    #: (data dir, kind) -> load_validation_map / load_diagnostic_map result.
    validations: Dict[Tuple[str, str], Dict] = field(default_factory=dict)
    diagnostics: Dict[Tuple[str, str], Dict] = field(default_factory=dict)

    def annotate(self, row: Dict[str, Any], data_dir: str, kind: str) -> None:
        """``Validation`` and ``Diagnostic`` (with badges) of a carrier row."""
        variant, sample = row.get("Variant", ""), row.get("Sample", "")
        add_validation_status_to_row(
            row, self.validations.get((data_dir, kind), {}), variant, sample
        )
        add_diagnostic_status_to_row(
            row, self.diagnostics.get((data_dir, kind), {}), variant, sample
        )
        row["DiagnosticBadges"] = diagnostic_badges(row)


def diagnostic_badges(row: Dict[str, Any]) -> List[Dict[str, str]]:
    """The row's diagnostic as a plain badge, curators in the tooltip.

    Not the tables' ``diagnostic`` cell: that one opens the diagnostic dialog
    on a click, and these pages only read. Diagnostics are recorded on the
    family pages.
    """
    data = row.get("Diagnostic_badge")
    if not data:
        return []
    badge = data["badge"]
    lines = [
        " ".join(p for p in (d["s"], "-", d["u"], d["t"], d["c"]) if p)
        for d in data["details"]
    ]
    return [
        {
            "label": row["Diagnostic"],
            "color": badge["bg"],
            "textColor": badge["tc"],
            "tooltip": "\n".join(lines),
        }
    ]


def load_curation(data_dirs: Iterable[str]) -> Curation:
    """Read every data directory's SNV, SV and STR curation files (blocking)."""
    cur = Curation()
    for d in data_dirs:
        for kind, vtype in _CURATED.items():
            cur.validations[(d, kind)] = load_validation_map(
                validation_path(Path(d), vtype), None
            )
            cur.diagnostics[(d, kind)] = load_diagnostic_map(
                diagnostic_path(Path(d), vtype), None
            )
    return cur


# ---------------------------------------------------------------------------
# Shared row pieces
# ---------------------------------------------------------------------------


def _carriers(df: pl.DataFrame, keys: List[str], joined: List[str]) -> pl.DataFrame:
    """One row per ``keys``: ``joined`` columns as their sorted distinct values,
    the others as their first value."""
    others = [c for c in df.columns if c not in keys and c not in joined]
    return (
        df.group_by(keys, maintain_order=True)
        .agg(
            *[
                pl.col(c).drop_nulls().unique().sort().str.join(", ").alias(c)
                for c in joined
            ],
            *[pl.col(c).first() for c in others],
        )
        .sort(keys)
    )


def _person(row: Dict[str, Any], dirs: Dict[str, DirInfo]) -> DirInfo | None:
    info = dirs.get(row.get("data_dir", ""))
    row["Data dir"] = info.label if info else row.get("data_dir", "")
    row["Cohort"] = row.get("cohort") or ""
    row["FID"] = row.get("family_id") or ""
    row["Sample"] = row.get("sample") or ""
    return info


def seen_summary(rows: List[Dict[str, Any]]) -> Dict[str, int]:
    """Individuals, families, cohorts and data directories among carrier rows."""
    return {
        "individuals": len({(r["data_dir"], r["Sample"]) for r in rows}),
        "families": len({(r["data_dir"], r["Cohort"], r["FID"]) for r in rows}),
        "cohorts": len({(r["data_dir"], r["Cohort"]) for r in rows}),
        "data_dirs": len({r["data_dir"] for r in rows}),
    }


def _gene_badges(value: Any, scorer: Any) -> List[Dict[str, Any]]:
    return build_gene_badges(split_genes(value), scorer)


# ---------------------------------------------------------------------------
# SNVs
# ---------------------------------------------------------------------------


def snv_carrier_rows(
    df: pl.DataFrame, curation: Curation, dirs: Dict[str, DirInfo], scorer: Any
) -> List[Dict[str, Any]]:
    """One row per carrier of each variant: data dir, cohort, family, sample.

    A carrier listed by several Wombat configs is one row naming them all.
    """
    if df.is_empty():
        return []
    frame = _carriers(
        df, ["data_dir", "cohort", "family_id", "sample", "variant"], ["config", "gene"]
    )
    rows = frame.to_dicts()
    unknown: Set[str] = set()
    for row in rows:
        info = _person(row, dirs)
        row["Variant"] = row["variant"]
        row["Genes"] = row["gene"]
        row["GeneBadges"] = _gene_badges(row["gene"], scorer)
        row["Configs"] = row["config"]
        row["VEP_Consequence"] = row.get("consequence") or ""
        row["VEP_CLIN_SIG"] = row.get("clin_sig") or ""
        add_vep_badges(row, unknown, unknown)
        row["gnomAD"] = row.get("gnomad")
        row["GT"] = row.get("sample_gt") or ""
        row["Father GT"] = row.get("father_gt") or ""
        row["Mother GT"] = row.get("mother_gt") or ""
        curation.annotate(row, row["data_dir"], "snv")
        first_config = (row["config"] or "").split(", ")[0]
        goto_fields(
            row, info, family_path(row["Cohort"], row["FID"], "snv", first_config)
        )
    return rows


def snv_variant_rows(carriers: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """One row per variant: how often, where, and what curators said of it."""
    by_variant: Dict[str, List[Dict[str, Any]]] = {}
    for row in carriers:
        by_variant.setdefault(row["Variant"], []).append(row)
    out: List[Dict[str, Any]] = []
    for variant, rows in by_variant.items():
        first = rows[0]
        seen = seen_summary(rows)
        statuses = Counter(r.get("Validation") or UNVALIDATED for r in rows)
        verdicts = sorted({r["Diagnostic"] for r in rows if r.get("Diagnostic")})
        out.append(
            {
                "Variant": variant,
                "Genes": first["Genes"],
                "GeneBadges": first["GeneBadges"],
                "ConsequenceBadges": first["ConsequenceBadges"],
                "ClinVarBadges": first["ClinVarBadges"],
                "VEP_Consequence": first["VEP_Consequence"],
                "VEP_CLIN_SIG": first["VEP_CLIN_SIG"],
                "gnomAD": first["gnomAD"],
                "Carriers": seen["individuals"],
                "Families": seen["families"],
                "Cohorts": ", ".join(sorted({r["Cohort"] for r in rows})),
                "Validations": validation_tally(statuses),
                "Diagnostics": ", ".join(verdicts),
                "_pos": _pos_key(variant),
            }
        )
    out.sort(key=lambda r: (-r["Carriers"], r["_pos"]))
    return out


def validation_tally(statuses: Counter) -> str:
    """``present 2 · TODO 1``, validated statuses first."""
    order = sorted(statuses, key=lambda s: (s == UNVALIDATED, s))
    return " · ".join(f"{s} {statuses[s]}" for s in order)


def _pos_key(variant: str) -> Tuple[str, int]:
    parts = variant.split(":")
    try:
        return parts[0], int(parts[1])
    except (IndexError, ValueError):
        return variant, 0


# ---------------------------------------------------------------------------
# SVs, STRs, ROH
# ---------------------------------------------------------------------------


def sv_rows(
    df: pl.DataFrame, curation: Curation, dirs: Dict[str, DirInfo], scorer: Any
) -> List[Dict[str, Any]]:
    """One row per call and carrier, the call's genes among the gene's names."""
    if df.is_empty():
        return []
    frame = _carriers(
        df, ["data_dir", "cohort", "family_id", "sample", "variant"], ["gene"]
    )
    rows = frame.to_dicts()
    for row in rows:
        info = _person(row, dirs)
        row["Variant"] = row["variant"]
        row["Locus"] = row["locus"]
        row["Type"] = row.get("sv_type") or ""
        row["wisecondorX"] = row.get("call") or ""
        try:
            row["Size (kb)"] = round((row["end"] - row["start"]) / 1000, 1)
        except TypeError:
            row["Size (kb)"] = None
        row["Ratio"] = row.get("ratio")
        row["Exonic"] = "yes" if row.get("exonic") else "no"
        row["GeneBadges"] = _gene_badges(row["gene"], scorer)
        curation.annotate(row, row["data_dir"], "sv")
        goto_fields(row, info, family_path(row["Cohort"], row["FID"], "sv"))
    return rows


def str_noteworthy(df: pl.DataFrame) -> pl.DataFrame:
    """Genotypes past normal: any pipeline severity, else a non-normal status.

    A table written before severities has none in any row and falls back to
    its statuses.
    """
    if df.is_empty():
        return df
    has_severity = pl.col("severity").is_not_null() & (pl.col("severity") != "")
    abnormal = pl.col("status").is_not_null() & (pl.col("status") != "normal")
    return df.filter(has_severity | (~has_severity & abnormal))


def str_rows(
    df: pl.DataFrame, curation: Curation, dirs: Dict[str, DirInfo], scorer: Any
) -> List[Dict[str, Any]]:
    """The noteworthy genotypes of the gene's loci, one row per member."""
    rows = str_noteworthy(df).to_dicts()
    for row in rows:
        info = _person(row, dirs)
        row["Variant"] = row["variant"]
        row["Locus"] = row["locus"]
        row["Catalog"] = row.get("label") or ""
        alleles = [a for a in (row.get("allele1"), row.get("allele2")) if a is not None]
        row["Alleles"] = "/".join(str(a) for a in alleles)
        status = status_info(row.get("status"))
        row["StatusBadges"] = [{"label": status["label"], "color": status["color"]}]
        severity = severity_info(row.get("severity"))
        row["SeverityBadges"] = (
            [{"label": severity["label"], "color": severity["color"]}]
            if severity["key"]
            else []
        )
        row["Disease"] = row.get("disease") or ""
        row["GeneBadges"] = _gene_badges(row["gene"], scorer)
        curation.annotate(row, row["data_dir"], "str")
        goto_fields(row, info, family_path(row["Cohort"], row["FID"], "str"))
    rows.sort(key=lambda r: (-severity_info(r.get("severity"))["rank"], r["Sample"]))
    return rows


def str_genotyped(df: pl.DataFrame) -> int:
    """Individuals genotyped at the gene's loci, noteworthy or not."""
    if df.is_empty():
        return 0
    return df.select(pl.struct("data_dir", "sample").n_unique()).item()


def roh_rows(df: pl.DataFrame, dirs: Dict[str, DirInfo]) -> List[Dict[str, Any]]:
    """One row per segment covering the gene."""
    if df.is_empty():
        return []
    frame = _carriers(
        df,
        ["data_dir", "cohort", "family_id", "sample", "chrom", "start", "end"],
        ["gene"],
    )
    rows = frame.to_dicts()
    for row in rows:
        info = _person(row, dirs)
        row["Segment"] = f"{row['chrom']}:{row['start']}-{row['end']}"
        row["Length (Mb)"] = row.get("length_mb")
        row["Class"] = row.get("roh_class") or ""
        fraction = row.get("gene_fraction")
        row["Gene covered"] = None if fraction is None else round(100 * fraction)
        goto_fields(row, info, family_path(row["Cohort"], row["FID"], "roh"))
    rows.sort(key=lambda r: -(r.get("length_mb") or 0))
    return rows


# ---------------------------------------------------------------------------
# Diagnostics and validations, straight from the curation files
# ---------------------------------------------------------------------------


def _read_curation_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def _family_cohorts(data_dir: str) -> Dict[str, str]:
    """Family -> its first cohort in the data directory."""
    store = get_all_data_stores().get(data_dir)
    out: Dict[str, str] = {}
    if store is not None:
        for name, cohort in sorted(store.cohorts.items()):
            for fid in cohort.families:
                out.setdefault(fid, name)
    return out


def diagnostic_rows(
    data_dirs: Iterable[str],
    dirs: Dict[str, DirInfo],
    scorer: Any,
    *,
    gene_keys: Optional[Set[str]] = None,
    variant: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """Standing diagnostics naming one of ``gene_keys`` (or of ``variant``).

    Read from the files, not the index: a diagnostic is found even when its
    variant is in no cohort table. One row per type, family, variant and
    sample, the curators' verdicts folded as on the Diagnostic pages.
    """
    out: List[Dict[str, Any]] = []
    for d in data_dirs:
        info = dirs.get(d)
        cohorts = _family_cohorts(d)
        for dtype, path in diagnostic_paths(Path(d)):
            groups: Dict[Tuple[str, str, str], List[Dict[str, str]]] = {}
            for raw in _read_curation_rows(path):
                if variant is not None and raw.get("Variant") != variant:
                    continue
                if gene_keys is not None:
                    # An SV's genes carry their overlap: "SCN2A:exonic"
                    genes = {
                        g.split(":", 1)[0].upper()
                        for g in split_genes(raw.get("Gene", ""), ",")
                    }
                    if not genes & gene_keys:
                        continue
                key = (
                    raw.get("FID", ""),
                    raw.get("Variant", ""),
                    raw.get("Sample", ""),
                )
                groups.setdefault(key, []).append(raw)
            for (fid, var, sample), raws in groups.items():
                row: Dict[str, Any] = {
                    "data_dir": d,
                    "Data dir": info.label if info else d,
                    "Type": dtype,
                    "FID": fid,
                    "Cohort": cohorts.get(fid, ""),
                    "Variant": var,
                    "Sample": sample,
                    "Genes": raws[-1].get("Gene", ""),
                    "GeneBadges": _gene_badges(raws[-1].get("Gene", ""), scorer),
                    "Impact": raws[-1].get("Impact", ""),
                    "ImpactBadges": impact_badges(dtype, raws[-1].get("Impact", "")),
                    "Users": ", ".join(sorted({r.get("User", "") for r in raws})),
                    "Timestamp": max(r.get("Timestamp", "") for r in raws),
                }
                diag_map = {
                    (var, sample): [
                        (
                            r.get("Diagnostic", ""),
                            r.get("User", ""),
                            r.get("Timestamp", ""),
                            r.get("Comment", ""),
                            r.get("Ignore", "0"),
                        )
                        for r in raws
                    ]
                }
                add_diagnostic_status_to_row(row, diag_map, var, sample)
                if not row["Diagnostic"]:  # every entry ignored
                    continue
                row["DiagnosticBadges"] = diagnostic_badges(row)
                goto_fields(row, info, f"/family/{fid}")
                out.append(row)
    out.sort(key=lambda r: r["Timestamp"], reverse=True)
    return out


def impact_badges(dtype: str, impact: str) -> List[Dict[str, str]]:
    """The diagnostic's ``Impact`` as the Diagnostic pages show it, roughly."""
    if not impact:
        return []
    if dtype == "SNV":
        return consequence_badges(impact, set())
    if dtype == "STR":
        info = status_info(impact)
        return [{"label": info["label"], "color": info["color"]}]
    return [{"label": impact, "color": "#6b7280"}]


def validation_history(
    data_dirs: Iterable[str], dirs: Dict[str, DirInfo], vtype: str, variant: str
) -> List[Dict[str, Any]]:
    """Every validation record of ``variant``, newest first, ignored ones marked.

    For the variant page: the full history, including a variant validated in a
    family whose results are in no cohort table.
    """
    out: List[Dict[str, Any]] = []
    for d in data_dirs:
        info = dirs.get(d)
        for raw in _read_curation_rows(validation_path(Path(d), vtype)):
            if raw.get("Variant") != variant:
                continue
            row = {
                "data_dir": d,
                "Data dir": info.label if info else d,
                "FID": raw.get("FID", ""),
                "Sample": raw.get("Sample", ""),
                "User": raw.get("User", ""),
                "Validation": raw.get("Validation", ""),
                "Inheritance": raw.get("Inheritance", ""),
                "Comment": raw.get("Comment", ""),
                "Ignored": "yes" if raw.get("Ignore") == "1" else "",
                "Timestamp": raw.get("Timestamp", ""),
            }
            goto_fields(row, info, f"/family/{row['FID']}")
            out.append(row)
    out.sort(key=lambda r: r["Timestamp"], reverse=True)
    return out
