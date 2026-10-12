"""Ancestry/PGS file discovery and the reference panel (utils/ancestry.py)."""

import pytest

from genetics_viz.utils import ancestry
from genetics_viz.utils.ancestry import (
    ROH_KINDS,
    _filename_pattern,
    _KINDS,
    _tag_strictness,
    _variant_key,
    find_ancestry_files,
    find_ancestry_files_in,
    find_cohort_ancestry_files,
    get_region_color,
    load_reference_pcs,
    probe_ancestry_data,
    probe_kind_in,
    probe_pgs_data,
)


@pytest.mark.parametrize(
    "name, expected",
    [
        ("F1.apgs_b1.0.0_dp10gq20.pcs.tsv", ("1.0.0", "dp10gq20", "pcs")),
        ("F1.apgs_b1.10.2_dp6gq15.pgs_zscore.tsv", ("1.10.2", "dp6gq15", "pgs_zscore")),
        ("F10.apgs_b1.0.0_dp10gq20.pcs.tsv", None),  # another family's prefix
        ("F1.apgs_b1.0.0_dp10gq20.pcs.tsv.gz", None),
        ("F1.apgs_b1.0.0.pcs.tsv", None),  # no filter tag
        ("F1.apgs_b1.0.0_dp10gq20.Q.tsv", None),  # kind not read
    ],
)
def test_filename_pattern(name, expected):
    m = _filename_pattern("F1", _KINDS).match(name)
    assert (m and (m["bundle"], m["filters"], m["kind"])) == expected


def test_roh_and_roh_genes_are_distinct_kinds():
    pat = _filename_pattern("F1", ROH_KINDS)
    assert pat.match("F1.apgs_b1.1.0_dp6gq15.roh.tsv")["kind"] == "roh"
    assert pat.match("F1.apgs_b1.1.0_dp6gq15.roh_genes.tsv")["kind"] == "roh_genes"
    assert pat.match("F1.apgs_b1.1.0_dp6gq15.froh.tsv")["kind"] == "froh"


def test_strictest_tag_ranks_first_not_string_order():
    tags = ["dp6gq15", "dp10gq20", "weird"]
    assert sorted(tags, key=_tag_strictness, reverse=True) == [
        "dp10gq20",
        "dp6gq15",
        "weird",
    ]


def test_newest_bundle_ranks_first():
    variants = [("1.9.0", "dp10gq20"), ("1.10.0", "dp6gq15"), ("1.10.0", "dp10gq20")]
    assert sorted(variants, key=_variant_key, reverse=True)[0] == ("1.10.0", "dp10gq20")
    assert sorted(variants, key=_variant_key, reverse=True)[-1] == ("1.9.0", "dp10gq20")


def test_release_ranks_above_its_release_candidate():
    assert _variant_key(("1.0.0", "dp10gq20")) > _variant_key(("1.0.0rc1", "dp10gq20"))
    assert _variant_key(("1.0.1rc1", "dp10gq20")) > _variant_key(("1.0.0", "dp10gq20"))


def test_find_files_selects_the_strictest_variant(data_dir):
    files = find_ancestry_files(data_dir, "1001")
    assert files.variant == ("1.1.0", "dp10gq20")
    assert files.all_variants == [("1.1.0", "dp10gq20"), ("1.1.0", "dp6gq15")]
    assert files.other_variants == [("1.1.0", "dp6gq15")]
    assert {"pcs", "ancestry", "pgs_zscore"} <= set(files.paths)
    assert files.label == "bundle 1.1.0 / dp10gq20"


def test_find_files_explicit_and_absent_variants(data_dir, tmp_path):
    assert (
        find_ancestry_files(data_dir, "1001", ("1.1.0", "dp6gq15")).filter_tag
        == "dp6gq15"
    )
    assert find_ancestry_files(data_dir, "1001", ("9.9.9", "dp6gq15")) is None
    assert find_ancestry_files_in(tmp_path / "nowhere", "X") is None


def test_cohort_files_and_probes(data_dir):
    class Cohort:
        name = "DEMO"
        path = data_dir / "cohorts" / "DEMO"

    assert find_cohort_ancestry_files(Cohort).pcs_path.name.startswith("DEMO.")
    assert probe_ancestry_data(data_dir, "1001") and probe_pgs_data(data_dir, "1001")
    assert not probe_ancestry_data(data_dir, "missing-family")
    assert not probe_kind_in(data_dir / "cohorts" / "DEMO" / "ancestry", "DEMO", "roh")


def test_region_colours_are_case_insensitive_with_a_fallback():
    assert get_region_color("eur") == get_region_color("EUR")
    assert (
        get_region_color(None)
        == get_region_color("Atlantis")
        == ancestry.REGION_COLORS["unknown"]
    )


def test_reference_panel_needs_configuration():
    assert load_reference_pcs("1.1.0") is None


def test_reference_panel_loads_and_is_cached(app_config, reference_root):
    ref = load_reference_pcs("1.1.0")
    assert ref.columns == ["IID", "PC1", "PC2", "PC3", "PC4", "population", "region"]
    assert len(ref) == 60
    assert load_reference_pcs("1.1.0") is ref


def test_missing_reference_bundle_is_none(app_config):
    assert load_reference_pcs("0.0.1") is None
