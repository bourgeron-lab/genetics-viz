"""Runs of homozygosity (ROH): file discovery, reference loading and readings.

The ghfc-ngs ``roh`` step runs ``ancestry-pgs roh`` on each family's ancestry
panel and writes, under ``<family>/roh/``, files named like the ancestry ones:
``<family_id>.apgs_b<bundle>_<filters>.<kind>.tsv`` for ``froh`` (one row per
individual), ``roh`` (segments), ``roh_genes`` and ``kinship``. The cohort
directory ``cohorts/<name>/roh/`` holds their concatenations with a
``family_id`` column, plus ``roh_windows`` - founders in ROH per 100 kb window.

Every individual is scored against unrelated reference genomes of their own
ancestry region (``ref_group``). Those reference values ship in the bundle's
``roh/`` tier, ``v<bundle>/bundle/roh/reference_froh.tsv.gz``, next to the
reference ROH hotspots. Plotting a person against that distribution is the
point of the display: the fixed published thresholds are crossed by a fifth
of unrelated South Asian reference genomes and half the Middle Eastern ones.

Thresholds, wording and colours come from ``config/roh.yaml``; nothing here or
in the tabs hard-codes them.
"""

import gzip
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import polars as pl
import yaml

from ..config_model import get_ancestry_reference_dir
from .ancestry import (
    NULL_VALUES,
    ROH_KINDS,
    AncestryFiles,
    Variant,
    find_ancestry_files_in,
    probe_kind_in,
)
from .cytobands import norm_chrom
from .sharding import get_family_path
from .tsv import read_tsv_or_none

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "roh.yaml"

#: Reference group meaning "every reference individual". The tool falls back
#: to it for outliers and low-confidence region labels.
POOLED = "ALL"

#: Derived column: population background, the sum of short and intermediate.
F_BACKGROUND = "f_short+intermediate"

#: Columns that hold identifiers or labels and must stay strings. Left to
#: inference, ``family_id`` and ``IID`` come out as integers when the first
#: rows happen to be numeric, and the read then fails on ``AU-ANG-FER-MPX3``.
_STRING_COLUMNS = {
    "IID",
    "IID1",
    "IID2",
    "family_id",
    "father",
    "mother",
    "chrom",
    "upd_candidate_chrom",
    "region",
    "af_source",
    "ref_group",
    "hotspot_group",
    "group",
    "gene",
    "gene_id",
    "roh_class",
    "froh_category",
    "homozygosity_pattern",
    "population",
}


# ---------------------------------------------------------------------------
# Display configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Grade:
    """One ``froh_category`` with its neutral grade and its relationship reading."""

    category: str
    grade: str
    color: str
    relationship: str
    threshold: Optional[float] = None


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


_CONFIG: Dict[str, Any] = _load_config()


def _grades() -> List[Grade]:
    return [
        Grade(
            category=g["category"],
            grade=g["grade"],
            color=g.get("color", "#9e9e9e"),
            relationship=g.get("relationship", ""),
            threshold=g.get("threshold"),
        )
        for g in _CONFIG.get("froh_grades", [])
    ]


GRADES: List[Grade] = _grades()
GRADE_HELP: str = _CONFIG.get("grade_help", "")
PATTERNS: Dict[str, Dict[str, str]] = _CONFIG.get("patterns", {})
CLASSES: Dict[str, Dict[str, str]] = _CONFIG.get("classes", {})
EXCESS_PERCENTILE: float = float(_CONFIG.get("excess_percentile", 95))
KINSHIP_DEGREES: List[Dict[str, Any]] = _CONFIG.get("kinship_degrees", [])


def reload_roh_config() -> None:
    """Reload the ROH display configuration from YAML."""
    global _CONFIG, GRADES, GRADE_HELP, PATTERNS, CLASSES, EXCESS_PERCENTILE
    global KINSHIP_DEGREES
    _CONFIG = _load_config()
    GRADES = _grades()
    GRADE_HELP = _CONFIG.get("grade_help", "")
    PATTERNS = _CONFIG.get("patterns", {})
    CLASSES = _CONFIG.get("classes", {})
    EXCESS_PERCENTILE = float(_CONFIG.get("excess_percentile", 95))
    KINSHIP_DEGREES = _CONFIG.get("kinship_degrees", [])


def grade_for(category: Optional[str]) -> Grade:
    """The grade of a ``froh_category``; an unknown one keeps its raw name."""
    for grade in GRADES:
        if grade.category == category:
            return grade
    name = category or "—"
    return Grade(category=name, grade=name, color="#9e9e9e", relationship=name)


def grade_text(category: Optional[str], show_relationship: bool) -> str:
    """Display text of a category: the grade, plus the reading for curators."""
    grade = grade_for(category)
    if show_relationship and grade.relationship and grade.category != "none":
        return f"{grade.grade} — {grade.relationship}"
    return grade.grade


def threshold_grades() -> List[Grade]:
    """Grades carrying a fixed threshold line, weakest first."""
    return sorted(
        (g for g in GRADES if g.threshold is not None), key=lambda g: g.threshold
    )


def pattern_info(pattern: Optional[str]) -> Dict[str, str]:
    """Label, colour and help text of a ``homozygosity_pattern``."""
    info = PATTERNS.get(pattern or "")
    if info:
        return info
    return {"label": pattern or "—", "color": "#9e9e9e", "help": ""}


def class_color(roh_class: Optional[str]) -> str:
    """Colour of a segment length class."""
    return CLASSES.get(roh_class or "", {}).get("color", "#9e9e9e")


def kinship_degree(phi: Optional[float]) -> Tuple[Optional[int], str]:
    """The KING degree of a kinship coefficient: ``(degree, label)``.

    The degree is ``None`` for unrelated, and the label "—" when ``phi`` is
    missing.
    """
    if phi is None:
        return None, "—"
    for entry in KINSHIP_DEGREES:
        if phi >= float(entry["min"]):
            return entry.get("degree"), entry["label"]
    return None, "Unrelated"


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------


def get_roh_dir(data_dir: Path, family_id: str) -> Path:
    """Return the ``roh/`` directory of a family (may not exist)."""
    return get_family_path(data_dir, family_id) / "roh"


def get_cohort_roh_dir(cohort: Any) -> Path:
    """Return the ``roh/`` directory of a cohort (may not exist)."""
    return cohort.path / "roh"


def find_roh_files(
    data_dir: Path, family_id: str, variant: Optional[Variant] = None
) -> Optional[AncestryFiles]:
    """Locate a family's ROH result files."""
    return find_ancestry_files_in(
        get_roh_dir(data_dir, family_id), family_id, variant, ROH_KINDS
    )


def find_cohort_roh_files(
    cohort: Any, variant: Optional[Variant] = None
) -> Optional[AncestryFiles]:
    """Locate a cohort's ROH result files."""
    return find_ancestry_files_in(
        get_cohort_roh_dir(cohort), cohort.name, variant, ROH_KINDS
    )


def probe_roh_data(data_dir: Path, family_id: str) -> bool:
    """Check whether a per-individual ROH table exists for this family."""
    return probe_kind_in(get_roh_dir(data_dir, family_id), family_id, "froh")


def probe_cohort_roh_data(cohort: Any) -> bool:
    """Check whether a cohort-level per-individual ROH table exists."""
    return probe_kind_in(get_cohort_roh_dir(cohort), cohort.name, "froh")


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def _header(path: Path) -> List[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def read_roh_table(path: Optional[Path]) -> Optional[pl.DataFrame]:
    """Read one ROH TSV with its identifier columns kept as strings.

    Returns ``None`` when the path is ``None`` or the file holds nothing. The
    ``chrom`` column is normalised (``chr1`` -> ``1``) so it joins with the
    reference hotspots and the cytobands, which carry no prefix.
    """
    if path is None:
        return None
    overrides = {c: pl.Utf8 for c in _header(path) if c in _STRING_COLUMNS}
    df = read_tsv_or_none(
        path,
        null_values=NULL_VALUES,
        infer_schema_length=10000,
        schema_overrides=overrides,
    )
    if df is None:
        return None
    if "chrom" in df.columns:
        df = df.with_columns(
            pl.col("chrom").map_elements(norm_chrom, return_dtype=pl.Utf8)
        )
    return df


# ---------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------

#: ``(file, bundle_version)`` -> frame, or None when absent. Both reference
#: tables are identical for every family, so they are read once per bundle.
_reference_cache: Dict[Tuple[str, str], Optional[pl.DataFrame]] = {}


def get_reference_roh_dir(bundle_version: str) -> Optional[Path]:
    """Return the bundle's ``roh/`` tier directory, or ``None`` when unconfigured.

    The returned path is not guaranteed to exist: bundles before 1.1.0 have no
    ROH tier, and the ROH step refuses them, so a missing directory beside real
    ROH results means the local reference copy is incomplete.
    """
    root = get_ancestry_reference_dir()
    if root is None:
        return None
    return root / f"v{bundle_version}" / "bundle" / "roh"


def _load_reference_file(bundle_version: str, name: str) -> Optional[pl.DataFrame]:
    directory = get_reference_roh_dir(bundle_version)
    if directory is None:
        return None
    path = directory / name
    key = (str(path), bundle_version)
    if key in _reference_cache:
        return _reference_cache[key]

    frame: Optional[pl.DataFrame] = None
    if not path.is_file():
        logger.warning(
            "ROH reference not found for bundle %s: %s", bundle_version, path
        )
    else:
        try:
            frame = read_roh_table(path)
            if frame is not None and len(frame) == 0:
                frame = None
        except Exception:
            logger.warning("Failed to read ROH reference %s", path, exc_info=True)
    _reference_cache[key] = frame
    return frame


def load_reference_froh(bundle_version: str) -> Optional[pl.DataFrame]:
    """ROH metrics of every unrelated reference individual, with their region.

    Columns: ``IID, population, region, af_source`` and the metric columns of
    ``froh.tsv`` (``froh_1.5mb``, ``f_short``, ...). ``None`` when the
    reference root is unset or the bundle has no ROH tier: a degraded state,
    the tabs then show the family's own values without the comparison.
    """
    return _load_reference_file(bundle_version, "reference_froh.tsv.gz")


def load_reference_hotspots(bundle_version: str) -> Optional[pl.DataFrame]:
    """ROH islands of the reference: ``group, chrom, start, end, max_freq``.

    ``group`` is a region or ``ALL``; ``chrom`` is normalised (no prefix).
    """
    return _load_reference_file(bundle_version, "hotspots.tsv.gz")


def clear_reference_cache() -> None:
    """Drop the cached reference tables."""
    _reference_cache.clear()


def reference_values(
    reference: pl.DataFrame, group: Optional[str], column: str
) -> List[float]:
    """The reference values of ``column`` for one comparison group.

    ``group`` follows the tool's ``ref_group``: a region name, or ``ALL`` for
    the pooled reference. :data:`F_BACKGROUND` is derived from its two parts.
    """
    frame = reference
    if group and group != POOLED:
        frame = frame.filter(pl.col("region") == group)
    if column == F_BACKGROUND:
        series = frame["f_short"] + frame["f_intermediate"]
    elif column in frame.columns:
        series = frame[column]
    else:
        return []
    return [float(v) for v in series.drop_nulls().to_list()]


def quantile(values: List[float], q: float) -> Optional[float]:
    """Linear-interpolated quantile of ``values`` (``q`` in 0-1)."""
    if not values:
        return None
    ordered = sorted(values)
    pos = q * (len(ordered) - 1)
    low = int(pos)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (pos - low)


def background(row: Dict[str, Any]) -> Optional[float]:
    """``f_short + f_intermediate`` of one ``froh`` row."""
    short, inter = to_float(row.get("f_short")), to_float(row.get("f_intermediate"))
    if short is None or inter is None:
        return None
    return short + inter


def to_float(value: Any) -> Optional[float]:
    """Coerce a cell to a float; ``None`` for missing, text or NaN."""
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return None if result != result else result


def metric_value(row: Dict[str, Any], column: str) -> Optional[float]:
    """A metric of one ``froh`` row, including the derived background."""
    if column == F_BACKGROUND:
        return background(row)
    return to_float(row.get(column))


# ---------------------------------------------------------------------------
# Kinship against the pedigree
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DeclaredRelation:
    """What the pedigree says about a pair."""

    label: str
    #: KING degree expected for the pair; ``None`` for a pair expected unrelated.
    degree: Optional[int]
    #: The two are the parents of someone in the family. Their kinship is the
    #: child's expected F, so it is a consanguinity reading, not a QC value.
    parents_pair: bool = False


def declared_relation(
    a: str, b: str, parents: Dict[str, Tuple[Optional[str], Optional[str]]]
) -> Optional[DeclaredRelation]:
    """The relationship the pedigree declares between ``a`` and ``b``, if any.

    ``parents`` maps an individual to ``(father, mother)``, either possibly
    ``None``.
    """
    fa, ma = parents.get(a, (None, None))
    fb, mb = parents.get(b, (None, None))
    if a in (fb, mb) or b in (fa, ma):
        return DeclaredRelation("Parent–offspring", 1)
    if fa and ma and fa == fb and ma == mb:
        return DeclaredRelation("Full sibs", 1)
    if (fa and fa == fb) or (ma and ma == mb):
        return DeclaredRelation("Half sibs", 2)
    for father, mother in parents.values():
        if {father, mother} == {a, b}:
            return DeclaredRelation("Parents of a child", None, parents_pair=True)
    return None


def kinship_rows(
    kinship: pl.DataFrame,
    parents: Dict[str, Tuple[Optional[str], Optional[str]]],
) -> List[Dict[str, Any]]:
    """Each kinship pair with its observed degree and the pedigree's declaration.

    ``mismatch`` is set when a pedigree-declared relative pair does not show
    the expected degree - a sample swap or a pedigree error. It is never set
    for the parents' pair or an undeclared pair: those carry no expectation.
    """
    rows: List[Dict[str, Any]] = []
    for row in kinship.to_dicts():
        a, b = row.get("IID1") or "", row.get("IID2") or ""
        phi = to_float(row.get("kinship"))
        degree, degree_label = kinship_degree(phi)
        declared = declared_relation(a, b, parents)
        mismatch = (
            declared is not None
            and not declared.parents_pair
            and declared.degree is not None
            and degree != declared.degree
        )
        rows.append(
            {
                "pair": f"{a} – {b}",
                "IID1": a,
                "IID2": b,
                "kinship": phi,
                "ibs0": to_float(row.get("ibs0")),
                "n_sites": row.get("n_sites"),
                "observed": degree_label,
                "declared": declared.label if declared else "Not declared",
                "is_declared": declared is not None and not declared.parents_pair,
                "parents_pair": bool(declared and declared.parents_pair),
                "mismatch": mismatch,
                "family_id": row.get("family_id"),
            }
        )
    return rows


def parents_from_froh(
    froh: pl.DataFrame,
) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
    """``IID -> (father, mother)`` from the pedigree columns ``froh`` carries."""
    out: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
    for row in froh.select(
        [c for c in ("IID", "father", "mother") if c in froh.columns]
    ).to_dicts():
        if row.get("IID"):
            out[row["IID"]] = (row.get("father"), row.get("mother"))
    return out


# ---------------------------------------------------------------------------
# Table badges, shared by the family summary and the cohort flagged list
# ---------------------------------------------------------------------------


def ordinal(value: Optional[float]) -> str:
    """``96.48`` -> ``"96th"``; "—" when missing."""
    if value is None:
        return "—"
    n = int(round(value))
    suffix = (
        "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    )
    return f"{n}{suffix}"


def grade_badge(category: Optional[str], show_relationship: bool) -> Dict[str, str]:
    """``badge_list`` entry for a ``froh_category``.

    Viewers get the neutral grade, with the grading rule as tooltip. Curators
    get the category name on the badge and the relationship reading on hover.
    """
    grade = grade_for(category)
    if show_relationship and grade.category != "none":
        return {
            "label": f"{grade.grade} · {grade.category}",
            "color": grade.color,
            "tooltip": grade.relationship,
        }
    return {"label": grade.grade, "color": grade.color, "tooltip": GRADE_HELP}


def pattern_badge(pattern: Optional[str]) -> Dict[str, str]:
    """``badge_list`` entry for a ``homozygosity_pattern``."""
    info = pattern_info(pattern)
    return {
        "label": info["label"],
        "color": info["color"],
        "tooltip": info.get("help", ""),
    }


def flag_badges(row: Dict[str, Any]) -> List[Dict[str, str]]:
    """Warning badges of one ``froh`` row: UPD candidate, failed QC."""
    badges = []
    upd = row.get("upd_candidate_chrom")
    if upd:
        badges.append(
            {
                "label": f"UPD? chr{upd}",
                "color": "#6a1b9a",
                "tooltip": (
                    "An isolated block of long ROH on one chromosome, with little "
                    "elsewhere: suggests uniparental isodisomy rather than parental "
                    "relatedness. This chromosome is left out of the grade."
                ),
            }
        )
    if row.get("qc_pass") is False:
        rate = to_float(row.get("call_rate"))
        badges.append(
            {
                "label": "Low call rate",
                "color": "#ef6c00",
                "tooltip": (
                    f"{rate:.1%} of panel sites called (< 90%). "
                    if rate is not None
                    else ""
                )
                + "Missing sites thin the evidence; short ROH become unreliable.",
            }
        )
    return badges


def is_flagged(row: Dict[str, Any]) -> bool:
    """True when a row deserves a look: graded, unusual pattern, UPD or failed QC."""
    return (
        (row.get("froh_category") or "none") != "none"
        or (row.get("homozygosity_pattern") or "typical") != "typical"
        or bool(row.get("upd_candidate_chrom"))
        or row.get("qc_pass") is False
    )
