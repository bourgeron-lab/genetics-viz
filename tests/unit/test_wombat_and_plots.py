"""Wombat aggregation and the pure plot builders (utils/wombat, ancestry_plots, roh_plots, ideogram)."""

import xml.etree.ElementTree as ET

import polars as pl

from genetics_viz.utils import ancestry_plots, ideogram, roh_plots
from genetics_viz.utils.ancestry import load_reference_pcs
from genetics_viz.utils.wombat import aggregate_by_variant


def wombat_rows():
    return pl.DataFrame(
        {
            "#CHROM": ["chr1", "chr1", "chr1", "chr2"],
            "POS": [10, 10, 10, 5],
            "REF": ["A", "A", "A", "G"],
            "ALT": ["T", "T", "T", "C"],
            "sample": ["S1", "S1", "S1", "S1"],
            "VEP_Consequence": [
                "missense_variant",
                "intron_variant",
                "missense_variant",
                ".",
            ],
            "VEP_SYMBOL": ["GENEB", "GENEA", "GENEB", None],
            "score": [0.5, None, 0.5, 1.0],
        }
    )


def test_aggregate_collapses_to_one_row_per_variant_and_sample():
    out = aggregate_by_variant(wombat_rows())
    assert out.height == 2 and out.columns[:6] == [
        "#CHROM",
        "POS",
        "REF",
        "ALT",
        "sample",
        "n_grouped",
    ]
    first = out.row(0, named=True)
    assert first["n_grouped"] == 3
    # Distinct values in first-seen order, so the string is stable.
    assert first["VEP_Consequence"] == "missense_variant,intron_variant"
    assert first["VEP_SYMBOL"] == "GENEB,GENEA"
    assert first["score"] == "0.5"
    second = out.row(1, named=True)
    assert second["VEP_Consequence"] == "" and second["VEP_SYMBOL"] == ""


def test_aggregate_is_deterministic():
    results = {tuple(aggregate_by_variant(wombat_rows()).row(0)) for _ in range(20)}
    assert len(results) == 1


# --- ancestry plots ----------------------------------------------------------------


def test_axis_bounds_contain_every_value():
    low, high = ancestry_plots.axis_bounds([0.0123, -0.0456, 0.031])
    assert low <= -0.0456 and high >= 0.031
    single = ancestry_plots.axis_bounds([0.5])
    assert single[0] < 0.5 < single[1]
    assert ancestry_plots.axis_bounds([]) is None
    assert ancestry_plots.round_coord("x") is None


def test_zoom_window_and_reference_series(app_config):
    ref = load_reference_pcs("1.1.0")
    rows = [{"PC1": 0.0, "PC2": 0.0}]
    window = ancestry_plots.zoom_window(rows, ref, ["CSA"], "PC1", "PC2")
    csa = ref.filter(pl.col("region") == "CSA")
    assert window[0][0] <= 0.0 and window[0][1] >= max(csa["PC1"])
    series = ancestry_plots.reference_series(ref, "PC1", "PC2")
    assert [s["name"] for s in series] == ["EUR", "CSA"]
    zoomed = ancestry_plots.reference_series(
        ref, "PC1", "PC2", window=((-1, 0.01), (-1, 1))
    )
    assert [s["name"] for s in zoomed] == ["EUR"]
    assert ancestry_plots.zoom_window([{"PC1": None}], None, [], "PC1", "PC2") is None


# --- ROH plots -----------------------------------------------------------------------


def test_axis_max_shows_the_next_threshold():
    assert roh_plots._axis_max([0.05], [0.03, 0.10, 0.17]) == 0.11
    assert roh_plots._axis_max([0.0], []) == 0.01
    assert roh_plots._axis_max([0.18], [0.03, 0.10, 0.17]) == 0.2


def test_reference_strip_is_deterministic_and_labelled():
    kwargs = dict(
        title="t",
        metric_label="m",
        p95=0.02,
        show_thresholds=True,
        show_relationship=False,
    )
    members = [{"name": "A", "value": 0.01}, {"name": "B", "value": 0.12}]
    first = roh_plots.reference_strip_option([0.01, 0.02, 0.03], members, **kwargs)
    assert first == roh_plots.reference_strip_option(
        [0.01, 0.02, 0.03], members, **kwargs
    )
    lines = first["series"][0]["markLine"]["data"]
    assert len(lines) == 4  # three thresholds and the 95th percentile
    assert not any(
        "degree" in line["name"] or "cousin" in line["name"] for line in lines
    )
    curator = roh_plots.reference_strip_option(
        [0.01], members, **{**kwargs, "show_relationship": True}
    )
    assert any(
        "cousin" in line["name"] for line in curator["series"][0]["markLine"]["data"]
    )


def test_genome_svg_is_well_formed_and_escaped():
    svg = roh_plots.genome_svg(
        ["<x>", "B"],
        {"<x>": [("2", 1_000_000, 9_000_000, "long", "seg <x>")]},
        {"B": [("1", 45_300_000, 46_200_000)]},
        {"<x>": "2"},
    )
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")
    assert "&lt;x&gt;" in svg and "<x>" not in svg
    assert 'stroke="#d32f2f"' in svg  # UPD outline


def test_islands_inserts_zeros_across_gaps_and_drops_non_autosomes():
    windows = [
        ("1", 0, 100_000, 0.5),
        ("1", 500_000, 600_000, 0.2),
        ("X", 0, 100_000, 0.9),
    ]
    option = roh_plots.islands_option(windows, [("1", 0, 100_000)], title="t")
    data = option["series"][1]["data"]
    xs = [p[0] for p in data]
    assert xs == sorted(xs)
    assert 0 in [p[1] for p in data[1:]]  # a dip to zero between the two windows
    assert all(p[2] != "X" for p in data)
    assert option["xAxis"]["max"] == round(
        sum(ideogram.chrom_length_mb(c) for c in roh_plots.AUTOSOMES), 1
    )


def test_chromosome_lengths_are_exact():
    assert 242.0 < ideogram.chrom_length_mb("2") < 242.3
    assert ideogram.chrom_length_mb("nope") == 0
