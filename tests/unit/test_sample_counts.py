"""Variants per sample and the cohort comparison (utils/sample_counts)."""

import os

import polars as pl
import pytest

import datagen
from genetics_viz.models import Cohort
from genetics_viz.utils import sample_counts
from genetics_viz.utils.sample_counts import (
    FOUNDERS,
    MissingHeaderError,
    NON_FOUNDERS,
    cohort_sample_counts,
    cohort_wcx_file,
    cohort_wombat_file,
    cohort_wombat_sources,
    compare_member,
    compare_to_cohort,
    count_variants_per_sample,
    family_header_sources,
    read_cohort_wombat,
    group_stats,
    histogram_bins,
)


def demo_cohort(data_dir):
    return Cohort.from_directory(data_dir / "cohorts" / "DEMO")


def test_cohort_files_are_found_by_exact_name(data_dir):
    path = cohort_wombat_file(data_dir, "DEMO", "vep1", "default")
    assert path is not None and path.name == "DEMO.rare.vep1.default.results.tsv"
    assert cohort_wombat_file(data_dir, "DEMO", "vep2", "default") is None
    assert cohort_wombat_file(data_dir, "DEMO", "vep1", "nope") is None
    # Older cohorts have only the name that keeps "annotated"
    strict = cohort_wombat_file(data_dir, "DEMO", "vep1", "strict")
    assert strict is not None and strict.name == (
        "DEMO.rare.vep1.annotated.strict.results.tsv"
    )
    assert cohort_wcx_file(data_dir, "DEMO").name == "DEMO_aberrations.bed"
    assert cohort_wcx_file(data_dir, "NOHEAD") is None


def test_counts_unique_variants_and_aberrations_per_sample(data_dir):
    wombat = cohort_wombat_file(data_dir, "DEMO", "vep1", "default")
    assert cohort_sample_counts(wombat, "wombat") == {"S3": len(datagen.WOMBAT_ROWS)}
    assert cohort_sample_counts(cohort_wcx_file(data_dir, "DEMO"), "wcx") == {
        "S3": 1,
        "S6": 1,
    }


def test_cohort_wombat_sources_are_sorted_and_listed_once(data_dir):
    wombat_dir = data_dir / "cohorts" / "DEMO" / "wombat"
    # Under both names, the table is listed once, as cohort_wombat_file finds it
    (wombat_dir / "DEMO.rare.vep1.annotated.default.results.tsv").write_text("")
    (wombat_dir / "DEMO.rare.vep0.default.results.tsv").write_text("")
    (wombat_dir / "OTHER.rare.vep1.default.results.tsv").write_text("")
    sources = cohort_wombat_sources(data_dir, "DEMO")
    assert [(s.wombat_config, s.vep_config, s.path.name) for s in sources] == [
        ("default", "vep0", "DEMO.rare.vep0.default.results.tsv"),
        ("default", "vep1", "DEMO.rare.vep1.annotated.default.results.tsv"),
        ("strict", "vep1", "DEMO.rare.vep1.annotated.strict.results.tsv"),
    ]
    assert cohort_wombat_sources(data_dir, "NOHEAD") == []


def test_a_headerless_cohort_table_is_read_with_a_family_header(data_dir, tmp_path):
    strict = cohort_wombat_file(data_dir, "DEMO", "vep1", "strict")
    sources = list(
        family_header_sources(data_dir, demo_cohort(data_dir), "vep1", "strict")
    )
    assert [p.name for p in sources] == ["1001.rare.vep1.annotated.strict.tsv"]
    with pytest.raises(MissingHeaderError, match="has no header line$"):
        read_cohort_wombat(strict)

    # A family whose file does not fit is passed over for one that does
    narrower = tmp_path / "narrower.tsv"
    narrower.write_text("\t".join(datagen.WOMBAT_COLUMNS[:-1]) + "\n")
    df = read_cohort_wombat(strict, [narrower, *sources], infer_schema_length=0)
    assert df.columns == datagen.WOMBAT_COLUMNS
    assert df["VEP_SYMBOL"].to_list() == ["SHANK3", "KIAA1109"]

    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    with pytest.raises(
        MissingHeaderError,
        match="14 columns against its 15; 1 other family file did not fit$",
    ):
        read_cohort_wombat(strict, [narrower, empty])


def test_a_headerless_cohort_file_borrows_the_family_header(data_dir, tmp_path):
    """Several real cohort de novo tables are a blank line, then rows."""
    strict = cohort_wombat_file(data_dir, "DEMO", "vep1", "strict")
    family = (
        datagen.family_dir(data_dir, "1001")
        / "wombat"
        / "1001.rare.vep1.annotated.strict.tsv"
    )
    # Counted without a header, every sample would come out at 0
    with pytest.raises(MissingHeaderError, match="has no header line"):
        cohort_sample_counts(strict, "wombat")
    assert sample_counts._count_cache == {}

    assert cohort_sample_counts(strict, "wombat", family) == {"S3": 2}

    narrower = tmp_path / "narrower.tsv"
    narrower.write_text("\t".join(datagen.WOMBAT_COLUMNS[:-1]) + "\n")
    sample_counts.clear_sample_count_cache()
    with pytest.raises(MissingHeaderError, match="columns against its"):
        cohort_sample_counts(strict, "wombat", narrower)


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
