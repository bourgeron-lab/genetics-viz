"""Search page logic without UI (utils/search_rows.py)."""

import polars as pl

from genetics_viz.utils import search_rows as sr
from genetics_viz.utils.view_presets import VIEW_PRESETS


def test_every_component_needs_a_locus_or_a_geneset():
    assert sr.missing_query([]) == "Please add at least one source component"
    ok = sr.new_wombat_component("wombat/x")
    ok["locus"]["value"] = "SHANK3"
    empty = sr.new_wcx_component("svs/wisecondorx")
    assert sr.missing_query([ok]) is None
    assert sr.missing_query([ok, empty]).startswith("Component 2 (wisecondorx)")


def test_wcx_calls_and_ratio_bounds():
    df = pl.DataFrame(
        {
            "wisecondorX": ["Robust LOSS", "Below threshold", "Robust GAIN"],
            "ratio": ["-1.2", "0.1", "0.4"],
        }
    )
    comp = sr.new_wcx_component()
    assert sr.filter_wcx_calls(df, comp)["wisecondorX"].to_list() == [
        "Robust LOSS",
        "Robust GAIN",
    ]
    comp["ratio_min"]["value"] = 0.5
    assert len(sr.filter_wcx_calls(df, comp)) == 1


def test_curated_coordinates_replace_the_locus():
    row = {"chr:start-end": "chr2:100-500", "_original_locus": "chr2:100-500"}
    vmap = {("k", "S1"): [("present", "", "", "0", "150", "", "t1")]}
    sr.apply_curated_coordinates(row, vmap, "k", "S1")
    assert row["Variant"] == "chr2:150-500" and row["IsCurated"]
    assert row["svlen"] == 350 and "Original: chr2:100-500" in row["_curated_tooltip"]
    # An ignored validation curates nothing.
    row = {"chr:start-end": "chr2:100-500"}
    sr.apply_curated_coordinates(
        row, {("k", "S1"): [("present", "", "", "1", "150", "", "t")]}, "k", "S1"
    )
    assert not row["IsCurated"] and row["chr:start-end"] == "chr2:100-500"


def test_individual_filters():
    ped = {
        "S1": {"Sex": "1", "Phenotype": "2", "Father": "S2", "Mother": "S3"},
        "S4": {"Sex": "2", "Phenotype": "1", "Father": "0", "Mother": ""},
    }
    rows = [{"sample": "S1"}, {"sample": "S4"}]
    f = sr.IndividualFilters()
    assert sr.filter_individuals(rows, ped, f) == rows
    f.has_parents["value"] = True
    assert sr.filter_individuals(rows, ped, f) == rows[:1]
    f.clear()
    f.exclude_samples["value"] = ["S1"]
    assert sr.filter_individuals(rows, ped, f) == rows[1:]


def test_sort_rows_puts_missing_values_last():
    rows = [{"x": "10"}, {"x": None}, {"x": "9"}]
    sr.sort_rows(rows, "x", "numerical", desc=False)
    assert [r["x"] for r in rows] == ["9", "10", None]


def test_validation_file_columns_follow_the_view():
    cols = ["Variant", "sample", "Validation", "ratio", "zscore", "type"]
    sv = sr.validation_file_columns("SV View", cols, VIEW_PRESETS)
    assert "Validation" not in sv and {"ratio", "zscore", "type"} <= set(sv)
    std = sr.validation_file_columns("Some View", cols, VIEW_PRESETS)
    assert "ratio" not in std or "ratio" in sr.preset_named(
        VIEW_PRESETS, "Standard View"
    ).get("columns", [])
