"""Gene sets for the filter menus: tiers, GeneTrek, local files."""

import polars as pl
from genetics_viz.utils import genesets
from genetics_viz.utils.gene_scoring import GeneScorer
from genetics_viz.utils.user_preferences import GenePrefs

COLORS = ("#c62828", "#ef6c00", "#f9a825", "#1565c0")


def test_sources_in_order(snapshot, data_dir):
    scorer = GeneScorer(snapshot, GenePrefs(COLORS, {"hc_ndd": 1, "constraint": 3}))
    sources = genesets.geneset_sources(data_dir, scorer)
    keys = [s.key for s in sources]
    assert keys[:2] == ["tier:1", "tier:3"]  # no tier 2 genes: not offered
    assert keys[-1] == "local:local_list"
    by_key = {s.key: s for s in sources}
    # The user's sets come before the rest of GeneTrek.
    assert keys.index("gt:hc_ndd") < keys.index("gt:syngo")
    assert by_key["gt:hc_ndd_v6"].secondary and not by_key["gt:hc_ndd"].secondary
    assert by_key["tier:1"].n_genes == 3
    assert "KIAA1109" in by_key["tier:1"].symbols  # old symbol expanded
    assert by_key["tier:3"].symbols >= {"SHANK3", "ARONLY"}


def test_selection_union(snapshot):
    scorer = GeneScorer(snapshot, GenePrefs(COLORS, {}))
    sources = genesets.geneset_sources(None, scorer)
    wanted = genesets.selected_symbols(sources, ["gt:syngo", "gt:constraint"])
    assert {"SHANK3", "ARONLY"} <= wanted and "BLTP1" not in wanted


def test_filter_expr_and_row_twin():
    df = pl.DataFrame(
        {"g": ["SHANK3", "a&kiaa1109", "X,Y", None, "", "SHANK", " plain ,Z"]}
    )
    wanted = frozenset({"SHANK3", "KIAA1109", "PLAIN"})
    expected = [True, True, False, False, False, False, True]
    got = df.select(genesets.geneset_filter_expr("g", wanted)).to_series().to_list()
    assert got == expected
    assert [genesets.cell_matches(v, wanted) for v in df["g"]] == expected


def test_local_genesets(data_dir):
    assert genesets.load_genesets(data_dir) == {"local_list": {"ARONLY", "PLAIN"}}
