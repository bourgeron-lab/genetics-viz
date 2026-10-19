"""Variants per sample, and where a family's members sit in their cohort.

The family Stats dialogs count each member's variants, which says little on its
own: in a trio there is nothing to compare against. So each member is also
placed on the cohort's per-sample distribution for the same Wombat config (or
the cohort's WisecondorX calls), read from the cohort-level file.

Both sides of that comparison are **unfiltered** counts from the cohort file:
the family table's filters (LCR, gnomAD, impact, genesets...) would have to be
replayed over a much larger table to be like-for-like, so they are not.

Founders and non-founders are compared within their own group. A de novo
config has rows only for children, so pooling them with the parents' zeros
would flag any child with two calls as an outlier.
"""

from __future__ import annotations

import bisect
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Tuple

import polars as pl

from genetics_viz.models import Cohort
from genetics_viz.utils.tsv import read_tsv_or_none

#: Tukey's fence: a count beyond Q3 + k·IQR (or below Q1 - k·IQR) is flagged.
OUTLIER_IQR_FACTOR = 1.5

#: The two populations a member is compared within.
FOUNDERS = "founders"
NON_FOUNDERS = "non-founders"

CountKind = Literal["wombat", "wcx"]


# --------------------------------------------------------------------------
# Cohort-level files
# --------------------------------------------------------------------------


def cohort_wombat_file(
    data_dir: Path, cohort_name: str, vep_config: str, wombat_config: str
) -> Optional[Path]:
    """The cohort's Wombat results for one VEP and Wombat config, if present.

    The family file is ``<FID>.rare.<vep>.annotated.<config>.tsv``. The
    cohort's adds ``results`` and, depending on the pipeline version, keeps
    ``annotated`` or drops it; older cohorts have only the annotated name.
    """
    wombat_dir = data_dir / "cohorts" / cohort_name / "wombat"
    for middle in (f"{vep_config}.annotated", vep_config):
        path = wombat_dir / f"{cohort_name}.rare.{middle}.{wombat_config}.results.tsv"
        if path.is_file():
            return path
    return None


def cohort_wcx_file(data_dir: Path, cohort_name: str) -> Optional[Path]:
    """The cohort's WisecondorX aberrations, plain or annotated, if present."""
    wcx_dir = data_dir / "cohorts" / cohort_name / "svs" / "wisecondorx"
    for name in (
        f"{cohort_name}_aberrations.bed",
        f"{cohort_name}_aberrations.annotated.bed",
    ):
        path = wcx_dir / name
        if path.is_file():
            return path
    return None


# --------------------------------------------------------------------------
# Counting
# --------------------------------------------------------------------------


def count_variants_per_sample(
    tsv_path: Path,
    geneset_genes: Set[str] | None = None,
) -> Dict[str, int]:
    """Count unique variants per sample in a cohort Wombat file.

    If geneset_genes is provided, only count variants where VEP_SYMBOL
    overlaps the gene set. Returns {sample_id: variant_count}.
    """
    df = read_tsv_or_none(
        tsv_path,
        infer_schema_length=0,
        null_values=[".", ""],
    )
    if df is None:
        return {}

    # Ensure required columns exist
    required = {"#CHROM", "POS", "REF", "ALT", "sample"}
    if not required.issubset(set(df.columns)):
        return {}

    # Filter by geneset genes if provided
    if geneset_genes and "VEP_SYMBOL" in df.columns:
        df = df.filter(
            pl.col("VEP_SYMBOL")
            .cast(pl.Utf8)
            .fill_null("")
            .str.split(",")
            .list.eval(pl.element().str.to_uppercase().is_in(list(geneset_genes)))
            .list.any()
        )

    # Count unique variants per sample
    variant_counts = (
        df.select(["#CHROM", "POS", "REF", "ALT", "sample"])
        .unique()
        .group_by("sample")
        .agg(pl.len().alias("count"))
    )

    return {
        str(row["sample"]): row["count"] for row in variant_counts.iter_rows(named=True)
    }


def count_aberrations_per_sample(bed_path: Path) -> Dict[str, int]:
    """Count WisecondorX aberrations per sample (the ``barcode`` column).

    Only the one column is needed, so the file is split by hand rather than
    through ``parse_wisecondorx_bed``, which classifies every call row by row.
    """
    counts: Dict[str, int] = {}
    with open(bed_path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        column = next((c for c in ("barcode", "sample") if c in header), None)
        if column is None:
            return {}
        idx = header.index(column)
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if idx < len(parts) and parts[idx]:
                counts[parts[idx]] = counts.get(parts[idx], 0) + 1
    return counts


_count_cache: Dict[Tuple[str, str, float, int], Dict[str, int]] = {}


def cohort_sample_counts(path: Path, kind: CountKind) -> Dict[str, int]:
    """Per-sample counts of a cohort file, cached on its (path, mtime, size).

    Blocking: call it through ``run.io_bound`` from a page.
    """
    stat = path.stat()
    key = (str(path), kind, stat.st_mtime, stat.st_size)
    cached = _count_cache.get(key)
    if cached is not None:
        return cached
    if kind == "wombat":
        counts = count_variants_per_sample(path)
    else:
        counts = count_aberrations_per_sample(path)
    _count_cache[key] = counts
    return counts


def clear_sample_count_cache() -> None:
    """Drop every cached per-sample count."""
    _count_cache.clear()


# --------------------------------------------------------------------------
# Comparison
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class GroupStats:
    """The per-sample count distribution of one cohort group."""

    name: str
    counts: List[int]  # sorted, one per sample, zero-filled
    median: float
    q1: float
    q3: float

    @property
    def n(self) -> int:
        return len(self.counts)

    @property
    def high_fence(self) -> float:
        return self.q3 + OUTLIER_IQR_FACTOR * (self.q3 - self.q1)

    @property
    def low_fence(self) -> float:
        return self.q1 - OUTLIER_IQR_FACTOR * (self.q3 - self.q1)


@dataclass(frozen=True)
class MemberComparison:
    """One family member's count against its cohort group."""

    sample_id: str
    count: int
    group: GroupStats
    #: Share of the group at or below this count, 0-100.
    percentile: float
    flag: Optional[Literal["high", "low"]]


def _quartiles(sorted_counts: List[int]) -> Tuple[float, float, float]:
    if len(sorted_counts) < 2:
        v = float(sorted_counts[0]) if sorted_counts else 0.0
        return v, v, v
    q1, median, q3 = statistics.quantiles(sorted_counts, n=4, method="inclusive")
    return q1, median, q3


def group_stats(name: str, counts: List[int]) -> GroupStats:
    """Median and quartiles of a group's per-sample counts."""
    ordered = sorted(counts)
    q1, median, q3 = _quartiles(ordered)
    return GroupStats(name=name, counts=ordered, median=median, q1=q1, q3=q3)


def compare_member(sample_id: str, count: int, group: GroupStats) -> MemberComparison:
    """Place one count on a group's distribution and flag it if out of range."""
    at_or_below = bisect.bisect_right(group.counts, count)
    percentile = 100.0 * at_or_below / group.n if group.n else 0.0
    flag: Optional[Literal["high", "low"]] = None
    if count > group.high_fence and count > group.median:
        flag = "high"
    elif count < group.low_fence and count < group.median:
        flag = "low"
    return MemberComparison(
        sample_id=sample_id,
        count=count,
        group=group,
        percentile=percentile,
        flag=flag,
    )


def compare_to_cohort(
    counts: Dict[str, int], cohort: Cohort, member_ids: List[str]
) -> Tuple[List[MemberComparison], List[str]]:
    """Compare family members with their cohort group.

    Every sample of the cohort pedigree is in the reference population, at 0
    when the file has no row for it. A sample in two families counts once.
    Returns the comparisons and the members absent from the cohort pedigree.
    """
    founder_of: Dict[str, bool] = {}
    for family in cohort.families.values():
        for sample in family.samples:
            # A sample in two families is a non-founder if it is one anywhere
            founder_of[sample.sample_id] = (
                founder_of.get(sample.sample_id, True) and sample.is_founder
            )

    by_group: Dict[str, List[int]] = {FOUNDERS: [], NON_FOUNDERS: []}
    for sid, is_founder in founder_of.items():
        by_group[FOUNDERS if is_founder else NON_FOUNDERS].append(counts.get(sid, 0))
    groups = {name: group_stats(name, vals) for name, vals in by_group.items()}

    comparisons: List[MemberComparison] = []
    missing: List[str] = []
    for sid in member_ids:
        if sid not in founder_of:
            missing.append(sid)
            continue
        group = groups[FOUNDERS if founder_of[sid] else NON_FOUNDERS]
        comparisons.append(compare_member(sid, counts.get(sid, 0), group))
    return comparisons, missing


def histogram_bins(counts: List[int], max_bins: int = 40) -> Tuple[List[int], int]:
    """Bin starts and bin width for a histogram of non-negative counts.

    Width 1 while the range fits in ``max_bins``, wider otherwise.
    """
    top = max(counts) if counts else 0
    width = max(1, -(-(top + 1) // max_bins))
    return list(range(0, top + 1, width)), width
