"""The denovo step on the family SNVs/Indels tab and the cohort De novo tab."""

import asyncio

import datagen
from nicegui import ui
from ui_helpers import badge_text, capture_tables, login


def _explain_not_parent(data_dir):
    """1001's NOT_PARENT explained: the de novo table is not held back."""
    with open(data_dir / "validations" / "sampleqc.tsv", "a") as f:
        f.write(
            "\t".join(["NOT_PARENT", "1001", "S1", "S3", "explained", "Adopted",
                       "curator", "2026-11-10T10:00:00"]) + "\n"
        )  # fmt: skip


def _variants(table):
    return [r["Variant"] for r in table["rows"]]


async def _rows_become(tables, variants, tries=100):
    for _ in range(tries):
        if tables and _variants(tables[-1]) == variants:
            return True
        await asyncio.sleep(0.05)
    return False


async def test_denovo_subtab_shows_strict_and_hidden_curation(
    user, data_dir, monkeypatch
):
    datagen.add_denovo_calls(data_dir)
    _explain_not_parent(data_dir)
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    # The De novo sub-tab comes first, from the step's Wombat config
    await user.should_see("Wombat Configuration: rare_denovo", retries=30)
    await user.should_see("De novo (1 strict / 2)", retries=10)
    await user.should_see(marker="denovo-children", retries=10)
    await user.should_see(marker="denovo-filter")
    # Strict: the flagged PLAIN call is hidden, and it was validated de novo
    assert await _rows_become(tables, ["chr22:50674000:G:A"])
    assert "dnm_flags" in tables[-1]["columns"]
    await user.should_see(marker="denovo-curation-banner", retries=10)
    await user.should_see("1 curated de novo call(s) hidden by the current filters")
    user.find(marker="denovo-show-hidden").click()
    assert await _rows_become(tables, ["chr22:50674000:G:A", "chr2:2000:G:T"])
    await user.should_not_see(marker="denovo-curation-banner", retries=10)
    # The cohort's verdict on S3 in the header
    badge = user.find(marker="denovo-qc-badge").elements.pop()
    assert badge.text == "De novo QC: S3"


async def test_flag_filter_hides_and_keeps(user, data_dir, monkeypatch):
    datagen.add_denovo_calls(data_dir, validation=False)
    _explain_not_parent(data_dir)
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001?tab=snvs&wcfg=denovo&wdm=all")
    assert await _rows_become(tables, ["chr22:50674000:G:A", "chr2:2000:G:T"])
    user.find(marker="denovo-flag-TR_HOMOPOLYMER").elements.pop().set_value("only")
    assert await _rows_become(tables, ["chr2:2000:G:T"])
    user.find(marker="denovo-flag-TR_HOMOPOLYMER").elements.pop().set_value("hide")
    assert await _rows_become(tables, ["chr22:50674000:G:A"])


async def test_old_links_open_the_renamed_tab(user, data_dir, monkeypatch):
    datagen.add_denovo_calls(data_dir, validation=False)
    _explain_not_parent(data_dir)
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001?tab=wombat&wcfg=rare_denovo")
    await user.should_see("Wombat Configuration: rare_denovo", retries=30)
    await user.should_see(marker="denovo-filter", retries=10)
    # The step's config has no sub-tab of its own beside De novo
    names = [t.props.get("name") for t in user.find(ui.tab).elements]
    assert "denovo" in names and "rare_denovo" not in names


async def test_without_the_step_the_subtab_says_so(user, data_dir):
    _explain_not_parent(data_dir)
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: rare_denovo", retries=30)
    await user.should_see(marker="denovo-not-run", retries=10)
    await user.should_not_see(marker="denovo-filter")
    await user.should_not_see(marker="denovo-qc-badge")


async def test_review_badge_counts_strict_candidates(user, data_dir):
    datagen.add_denovo_calls(data_dir, validation=False)
    _explain_not_parent(data_dir)
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    assert await badge_text(user, "snvs-tab-count") == "1"


async def test_variant_dialog_shows_the_evidence(user, data_dir, monkeypatch):
    datagen.add_denovo_calls(data_dir, validation=False)
    _explain_not_parent(data_dir)
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    assert await _rows_become(tables, ["chr22:50674000:G:A"])
    table = tables[-1]
    row = dict(table["rows"][0])
    # What the browser sends back: the row's own cells, not the hidden ones
    sent = {k: v for k, v in row.items() if not k.startswith("_")}
    with user.client:
        table["on_row_action"]({"action": "view_variant", "row": sent})
    await user.should_see(marker="denovo-evidence", retries=20)
    await user.should_see("De novo evidence")


async def test_cohort_tab_compares_with_the_literature(user, data_dir, monkeypatch):
    datagen.add_denovo_calls(data_dir)
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO?tab=denovo")
    await user.should_see(marker="denovo-references", retries=30)
    await user.should_see("Implied SNV mutation rate (per bp per generation)")
    await user.should_see("Jónsson 2017 (80.4%)", retries=5)
    children = next(t for t in tables if "n_strict" in t["columns"])
    assert [(r["child"], r["qc"]) for r in children["rows"]] == [("S3", "low")]
    flags = next(t for t in tables if "alone" in t["columns"])
    assert {r["flag"] for r in flags["rows"]} == {
        "TR_HOMOPOLYMER",
        "PARENT_ALT_POSSIBLE",
    }
    issues = next(t for t in tables if "problem_label" in t["columns"])
    assert [(r["Variant"], r["problem"]) for r in issues["rows"]] == [
        ("chr2:2000:G:T", "hidden")
    ]
    charts = [c.options["title"]["text"] for c in user.find(ui.echart).elements]
    assert "Strict candidates per callable Gb" in charts
    # S3: QC low (medium colour), and the hidden curated call
    for _ in range(100):
        try:
            pills = [e.text for e in user.find(marker="denovo-tab-pills").elements]
        except AssertionError:
            pills = []
        if pills:
            break
        await asyncio.sleep(0.05)
    assert sorted(pills) == ["1", "1"]


async def test_cohort_tab_reader_sees_no_curation(user, data_dir, monkeypatch):
    datagen.add_denovo_calls(data_dir)
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO?tab=denovo")
    await user.should_see(marker="denovo-references", retries=30)
    assert not any("problem_label" in t["columns"] for t in tables)
    flags = next(t for t in tables if "alone" in t["columns"])
    assert "present" not in flags["columns"]


async def test_cohort_tab_disabled_without_tables(user, data_dir):
    await login(user, "curator")
    await user.open("/cohort/DEMO")
    await user.should_see("Individuals", retries=10)
    tab = next(
        t for t in user.find(ui.tab).elements if t.props.get("name") == "De novo"
    )
    assert tab.props.get("disable") is True


async def test_cell_line_child_gets_the_banner(user, data_dir, monkeypatch):
    # S3 a cell line, its count high and not germline-like (genome QC).
    datagen.add_denovo_calls(data_dir, validation=False, genome_qc=True)
    _explain_not_parent(data_dir)
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: rare_denovo", retries=30)
    await user.should_see(marker="denovo-germline-banner", retries=10)
    await user.should_see("S3: DNA from a cell line (expected germline share 4%")
    await user.should_see("Cell line", retries=5)
    badge = user.find(marker="denovo-qc-badge").elements.pop()
    assert badge.text == "De novo QC: S3"


async def test_cohort_tab_leaves_cell_lines_out_of_the_literature(
    user, data_dir, monkeypatch
):
    datagen.add_denovo_calls(data_dir, validation=False, genome_qc=True)
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO?tab=denovo")
    await user.should_see(marker="denovo-references", retries=30)
    children = next(t for t in tables if "n_strict" in t["columns"])
    (row,) = children["rows"]
    assert row["cell_line"] and row["qc_flags"] == "CELL_LINE,NOT_GERMLINE_LIKE"
    assert {"qc_flags", "seq_class", "expected_germline_share"} <= set(
        children["columns"]
    )
    # The only child is a cell line: nothing left to compare with the
    # literature, and the note says why.
    note = user.find(marker="denovo-cohort-note").elements.pop().text
    assert "without the 1 cell-line child (S3)" in note
    assert "0 strict on the autosomes" in note
    chart = next(
        c
        for c in user.find(ui.echart).elements
        if c.options["title"]["text"] == "Strict candidates per callable Gb"
    )
    assert chart.options["series"][0]["data"][0]["itemStyle"]["color"] == "#880e4f"
    assert "purple: cell lines" in chart.options["title"]["subtext"]
