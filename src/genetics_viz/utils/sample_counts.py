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
import re
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Literal, Optional, Set, Tuple

import polars as pl

from genetics_viz.models import Cohort
from genetics_viz.utils.sharding import get_family_path
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


@dataclass(frozen=True)
class CohortWombatSource:
    """One cohort Wombat results table, as the Search page offers it."""

    path: Path
    vep_config: str
    wombat_config: str


def cohort_wombat_sources(data_dir: Path, cohort_name: str) -> List[CohortWombatSource]:
    """Every Wombat results table of a cohort, sorted by config, then VEP config.

    Sorted so the first, the Search page's default, does not depend on the
    order the directory lists in. A table under both names (with and without
    ``annotated``) is listed once, under the name ``cohort_wombat_file`` picks.
    """
    wombat_dir = data_dir / "cohorts" / cohort_name / "wombat"
    if not wombat_dir.is_dir():
        return []
    pattern = re.compile(
        rf"{re.escape(cohort_name)}\.rare\.([^.]+)\.(?:annotated\.)?(.+?)\.results\.tsv$"
    )
    configs = {
        (m.group(2), m.group(1))
        for f in wombat_dir.glob("*.results.tsv")
        if (m := pattern.match(f.name))
    }
    sources: List[CohortWombatSource] = []
    for wombat_config, vep_config in sorted(configs):
        path = cohort_wombat_file(data_dir, cohort_name, vep_config, wombat_config)
        if path is not None:
            sources.append(CohortWombatSource(path, vep_config, wombat_config))
    return sources


def family_wombat_file(
    data_dir: Path, family_id: str, vep_config: str, wombat_config: str
) -> Path:
    """Where a family's Wombat table for one VEP and Wombat config would be."""
    return (
        get_family_path(data_dir, family_id)
        / "wombat"
        / f"{family_id}.rare.{vep_config}.annotated.{wombat_config}.tsv"
    )


def family_header_sources(
    data_dir: Path, cohort: Cohort, vep_config: str, wombat_config: str
) -> Iterator[Path]:
    """The cohort's family files of one config, to lend a headerless table theirs.

    Lazy, in family order: a header that fits is usually found in the first
    family, and each existence check is a round trip on the data mount.
    """
    for family_id in sorted(cohort.families):
        path = family_wombat_file(data_dir, family_id, vep_config, wombat_config)
        if path.is_file():
            yield path


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


class MissingHeaderError(ValueError):
    """A cohort Wombat file with no header line, and none that fits to borrow."""


def _first_lines(path: Path) -> Tuple[int, Optional[str]]:
    """The number of leading blank lines, and the first non-blank line."""
    blanks = 0
    with open(path, "r") as f:
        for line in f:
            if line.strip():
                return blanks, line.rstrip("\n")
            blanks += 1
    return blanks, None


def _wombat_read_kwargs(
    tsv_path: Path, header_sources: Iterable[Path]
) -> Optional[Dict[str, Any]]:
    """``read_csv`` arguments for a cohort Wombat file, or None when it is blank.

    Some cohort de novo files have no header: a blank line, then the rows. The
    first of ``header_sources`` (family files of the same config) whose header
    has as many columns as the rows lends it; with none, ``MissingHeaderError``.
    """
    blanks, first = _first_lines(tsv_path)
    if first is None:
        return None
    read_kwargs: Dict[str, Any] = {"skip_rows": blanks}
    if first.startswith("#CHROM"):
        return read_kwargs
    width = len(first.split("\t"))
    tried: List[Tuple[Path, int]] = []
    for source in header_sources:
        with open(source, "r") as f:
            borrowed = f.readline().rstrip("\n").split("\t")
        if len(borrowed) == width and borrowed[0] == "#CHROM":
            read_kwargs.update(has_header=False, new_columns=borrowed)
            return read_kwargs
        # A zero-byte family file (nothing to report) reads as one column
        tried.append((source, len(borrowed)))
    detail = ""
    if tried:
        source, columns = tried[0]
        detail = f", and {source.name} has {columns} columns against its {width}"
        others = len(tried) - 1
        if others:
            detail += f"; {others} other family file{'s' * (others > 1)} did not fit"
    raise MissingHeaderError(f"{tsv_path.name} has no header line{detail}")


def read_cohort_wombat(
    tsv_path: Path, header_sources: Iterable[Path] = (), **kwargs: Any
) -> Optional[pl.DataFrame]:
    """Read a cohort Wombat table, borrowing a family's header if it has none.

    ``header_sources`` are family files of the same config, tried in turn
    (see :func:`family_header_sources`); they are read only when the cohort
    file has no header. Raises ``MissingHeaderError`` when none fits, rather
    than reading the first row as column names. Other keyword arguments go to
    :func:`read_tsv_or_none`. None for a file with no rows.
    """
    read_kwargs = _wombat_read_kwargs(tsv_path, header_sources)
    if read_kwargs is None:
        return None
    return read_tsv_or_none(tsv_path, **kwargs, **read_kwargs)


def count_variants_per_sample(
    tsv_path: Path,
    geneset_genes: Set[str] | None = None,
    header_sources: Iterable[Path] = (),
) -> Dict[str, int]:
    """Count unique variants per sample in a cohort Wombat file.

    If geneset_genes is provided, only count variants where VEP_SYMBOL
    overlaps the gene set. Returns {sample_id: variant_count}.

    ``header_sources`` (family files of the same config, tried in turn) lend
    their header to a cohort file that has none; see :func:`read_cohort_wombat`.
    Without one that fits, ``MissingHeaderError`` is raised rather than
    counting every sample as 0.
    """
    df = read_cohort_wombat(
        tsv_path,
        header_sources,
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


def cohort_sample_counts(
    path: Path, kind: CountKind, header_source: Optional[Path] = None
) -> Dict[str, int]:
    """Per-sample counts of a cohort file, cached on its (path, mtime, size).

    ``header_source`` is a family Wombat file of the same config, for a cohort
    file written without a header. Its families share the columns, so the
    cache key leaves it out. Blocking: call it through ``run.io_bound``.
    """
    stat = path.stat()
    key = (str(path), kind, stat.st_mtime, stat.st_size)
    cached = _count_cache.get(key)
    if cached is not None:
        return cached
    if kind == "wombat":
        sources = [header_source] if header_source is not None else []
        counts = count_variants_per_sample(path, header_sources=sources)
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
