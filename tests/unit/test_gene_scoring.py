"""Gene badges: best tier wins, tooltips, links, caps."""

from genetics_viz.utils import gene_scoring as gs
from genetics_viz.utils.user_preferences import GenePrefs

COLORS = ("#c62828", "#ef6c00", "#f9a825", "#1565c0")


def scorer(snapshot, sets) -> gs.GeneScorer:
    return gs.GeneScorer(snapshot, GenePrefs(COLORS, sets))


def test_best_tier_wins(snapshot):
    s = scorer(snapshot, {"hc_ndd": 2, "constraint": 1, "syngo": 4})
    assert s.tier("SHANK3") == 1 and s.color("SHANK3") == COLORS[0]
    assert s.tier("BLTP1") == 2
    assert s.tier("KIAA1109") == 2  # old symbol, same gene
    assert s.rank("SHANK3") > s.rank("BLTP1") > s.rank("PLAIN") == 0


def test_show_and_hidden_do_not_colour(snapshot):
    s = scorer(snapshot, {"hc_ndd": "show"})
    assert s.tier("SHANK3") is None and s.color("SHANK3") == gs.NEUTRAL_COLOR


def test_versioned_and_stale_keys(snapshot):
    s = scorer(snapshot, {"hc_ndd_v6": 3, "gone_list": 1})
    assert s.tier("SHANK3") == 3
    assert s.stale_keys == ["gone_list"]


def test_text_colour_by_luminance():
    assert gs.text_color_for("#ffffff") == "black"
    assert gs.text_color_for("#f9a825") == "black"
    assert gs.text_color_for("#c62828") == "white"
    assert gs.text_color_for("#1565c0") == "white"
    assert gs.text_color_for("garbage") == "black"


def test_tooltip(snapshot):
    s = scorer(snapshot, {"hc_ndd": 1, "constraint": "show"})
    tip = s.tooltip("SHANK3")
    assert tip.splitlines()[0] == "SHANK3 — SH3 and ankyrin 3"
    assert "Tier 1: High-confidence NDD v7" in tip
    assert "Also in: Constrained genes" in tip
    assert "+ 3 other GeneTrek sets, click to open" in tip
    assert "EAGLE Definitive" in tip and "Inheritance:" in tip
    assert s.tooltip("KIAA1109").startswith("KIAA1109 → BLTP1")
    assert s.tooltip("NOPE") == "NOPE: not in GeneTrek"
    assert "In none of your genesets" in s.tooltip("PLAIN")


def test_tooltip_without_data():
    from genetics_viz.utils.genetrek import GenetrekSnapshot

    s = gs.GeneScorer(GenetrekSnapshot(error="boom"), GenePrefs(COLORS, {}))
    assert "unavailable: boom" in s.tooltip("SHANK3")
    assert s.color("SHANK3") == gs.NEUTRAL_COLOR


def test_badges_sorted_capped_linked(snapshot, app_config):
    s = scorer(snapshot, {"hc_ndd": 2, "constraint": 1})
    genes = ["PLAIN", "BLTP1", "SHANK3", "A", "B", "C", "D", "E"]
    badges = gs.build_gene_badges(genes, s, cap=6)
    assert [b["label"] for b in badges[:2]] == ["SHANK3", "BLTP1"]
    assert len(badges) == 7 and badges[-1]["label"] == "+2 genes"
    assert badges[-1]["href"].startswith("http://genetrek.test/#search=")
    assert badges[0]["href"] == "http://genetrek.test/?search=SHANK3"
    assert gs.badge_labels(badges) == [b["label"] for b in badges[:6]]


def test_badge_options(snapshot, app_config):
    s = scorer(snapshot, {"hc_ndd": 1})
    badges = gs.build_gene_badges(
        ["PLAIN", "SHANK3"], s, cap=None, sort=False, bordered={"SHANK3"}
    )
    assert [b["label"] for b in badges] == ["PLAIN", "SHANK3"]
    assert "borderColor" in badges[1] and "borderColor" not in badges[0]
    assert badges[0]["href"].endswith("?search=PLAIN")


def test_split_genes_and_html(snapshot, app_config):
    assert gs.split_genes("A&B,C, A ,") == ["A", "B", "C"]
    assert gs.split_genes(None) == []
    badge = gs.gene_badge("SHANK3", scorer(snapshot, {}))
    badge["tooltip"] = 'say "<hi>"'
    html = gs.gene_badge_html(badge)
    assert html.startswith('<a class="q-badge gene-badge" href="http://genetrek.test/')
    assert "&quot;&lt;hi&gt;&quot;" in html


def test_scorer_cache_per_user(snapshot, app_config):
    a = gs.get_gene_scorer("reader")
    assert gs.get_gene_scorer("reader") is a
    assert gs.get_gene_scorer("curator") is not a
