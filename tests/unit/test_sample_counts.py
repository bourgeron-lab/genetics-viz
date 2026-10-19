"""Variants per sample and the cohort comparison (utils/sample_counts)."""

import os

import polars as pl

import datagen
from genetics_viz.models import Cohort
from genetics_viz.utils import sample_counts
from genetics_viz.utils.sample_counts import (
    FOUNDERS,
    NON_FOUNDERS,
    cohort_sample_counts,
    cohort_wcx_file,
    cohort_wombat_file,
    compare_member,
    compare_to_cohort,
    count_variants_per_sample,
    group_stats,
    histogram_bins,
)


def demo_cohort(data_dir):
    return Cohort.from_directory(data_dir / "cohorts" / "DEMO")


def test_cohort_files_are_found_by_exact_name(data_dir):
    path = cohort_wombat_file(data_dir, "DEMO", "vep1", "default")
    assert path is not None and path.name == "DEMO.rare.vep1.default.results.tsv"
    assert cohort_wombat_file(data_dir, "DEMO", "vep2", "default") is None
    assert cohort_wombat_file(data_dir, "DEMO", "vep1", "strict") is None
    # Older cohorts have only the name that keeps "annotated"
    annotated = path.with_name("DEMO.rare.vep1.annotated.strict.results.tsv")
    annotated.write_text(path.read_text())
    assert cohort_wombat_file(data_dir, "DEMO", "vep1", "strict") == annotated
    assert cohort_wcx_file(data_dir, "DEMO").name == "DEMO_aberrations.bed"
    assert cohort_wcx_file(data_dir, "NOHEAD") is None


def test_counts_unique_variants_and_aberrations_per_sample(data_dir):
    wombat = cohort_wombat_file(data_dir, "DEMO", "vep1", "default")
    assert cohort_sample_counts(wombat, "wombat") == {"S3": len(datagen.WOMBAT_ROWS)}
    assert cohort_sample_counts(cohort_wcx_file(data_dir, "DEMO"), "wcx") == {
        "S3": 1,
        "S6": 1,
    }


def test_cache_picks_up_a_rewritten_file(data_dir):
    path = cohort_wcx_file(data_dir, "DEMO")
    assert cohort_sample_counts(path, "wcx")["S3"] == 1
    with open(path, "a") as f:
        f.write("\t".join(map(str, datagen.WCX_ROWS[0])) + "\n")
    stat = path.stat()
    os.utime(path, (stat.st_atime, stat.st_mtime + 10))
    assert cohort_sample_counts(path, "wcx")["S3"] == 2
    sample_counts.clear_sample_count_cache()
    assert sample_counts._count_cache == {}


def test_members_are_compared_within_founders_or_non_founders(data_dir):
    comparisons, missing = compare_to_cohort(
        {"S3": 6}, demo_cohort(data_dir), ["S1", "S3", "NOT_THERE"]
    )
    assert missing == ["NOT_THERE"]
    s1, s3 = comparisons
    # S1, S2, S4, S5 are founders, all zero-filled
    assert s1.group.name == FOUNDERS and s1.group.counts == [0, 0, 0, 0]
    assert s1.count == 0 and s1.flag is None and s1.percentile == 100
    # S3 and S6 are non-founders (S6 has one parent in the pedigree)
    assert s3.group.name == NON_FOUNDERS and s3.group.counts == [0, 6]
    assert s3.percentile == 100


def test_flags_outliers_beyond_the_tukey_fence():
    group = group_stats(NON_FOUNDERS, [1, 2, 2, 2, 3, 3, 3, 40])
    assert group.median == 2.5 and group.q1 == 2 and group.q3 == 3
    assert compare_member("A", 40, group).flag == "high"
    assert compare_member("B", 3, group).flag is None
    assert compare_member("C", 0, group).flag == "low"


def test_with_no_spread_any_departure_from_the_bulk_is_flagged():
    group = group_stats(FOUNDERS, [0] * 9 + [1])
    assert group.q1 == group.q3 == 0
    assert compare_member("A", 1, group).flag == "high"
    assert compare_member("B", 0, group).flag is None


def test_histogram_bins_widen_past_the_maximum_bin_count():
    assert histogram_bins([0, 3]) == ([0, 1, 2, 3], 1)
    starts, width = histogram_bins([0, 99], max_bins=10)
    assert width == 10 and starts[-1] == 90
    assert histogram_bins([]) == ([0], 1)


def test_counts_unique_variants_in_a_geneset(tmp_path):
    path = tmp_path / "w.tsv"
    pl.DataFrame(
        {
            "#CHROM": ["1", "1", "2"],
            "POS": [1, 1, 2],
            "REF": ["A", "A", "C"],
            "ALT": ["T", "T", "G"],
            "sample": ["S1", "S1", "S2"],
            "VEP_SYMBOL": ["GENEA", "GENEA", "GENEB,GENEC"],
        }
    ).write_csv(path, separator="\t")
    assert count_variants_per_sample(path, None) == {"S1": 1, "S2": 1}
    assert count_variants_per_sample(path, {"GENEC"}) == {"S2": 1}
