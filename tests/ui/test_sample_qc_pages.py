"""The Sample QC views, in NiceGUI's user simulation, logged in as each role."""

from nicegui import ui
from ui_helpers import capture_tables, login, open_tab


def _tab(user, name):
    return next(t for t in user.find(ui.tab).elements if t.props.get("name") == name)


async def open_family_qc(user, monkeypatch, username, fid="1001"):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open(f"/cohort/DEMO/family/{fid}")
    await user.should_see(f"Family: {fid}", retries=10)
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    # The cohort's tables come in off-thread and add the cross-family flag.
    await user.should_not_see("Reading the cohort's tables", retries=30)
    return tables


def _issues(tables):
    return [t for t in tables if "severity" in t["columns"]][-1]


async def test_curator_sees_flags_evidence_and_cross_family_relatives(
    user, monkeypatch
):
    tables = await open_family_qc(user, monkeypatch, "curator")
    issues = _issues(tables)
    assert "evidence" in issues["columns"] and "other_family" in issues["columns"]
    assert [r["flag"] for r in issues["rows"]] == [
        "Not the declared parent",
        "Related across families",
        "Possible contamination",
    ]
    cross = next(r for r in issues["rows"] if r["other_family"])
    assert cross["other_family"] == "AU-X-1" and cross["samples"] == "S2 – S5"

    pairs = [t for t in tables if "expected_kinship" in t["columns"]][-1]
    assert "observed" in pairs["columns"]
    # Mismatches first: the declared father, then the relative in AU-X-1.
    assert [r["pair"] for r in pairs["rows"][:2]] == ["S1 – S3", "S2 – S5"]
    assert {r["pair"] for r in pairs["rows"]} >= {"S1 – S2"}  # the parents' pair

    sex = [t for t in tables if "karyotype" in t["columns"]][-1]
    assert {r["karyotype"] for r in sex["rows"]} == {"XX", "XY"}
    charts = [c.options for c in user.find(ui.echart).elements]
    assert any(s.get("type") == "heatmap" for c in charts for s in c.get("series", []))


async def test_reader_gets_neutral_wording_and_declared_pairs_only(user, monkeypatch):
    tables = await open_family_qc(user, monkeypatch, "reader")
    issues = _issues(tables)
    assert "evidence" not in issues["columns"]
    assert "Pedigree relationship inconsistent" in {r["flag"] for r in issues["rows"]}
    assert all(
        "S1 is not" not in b["tooltip"]
        for r in issues["rows"]
        for b in r["flag_badges"]
    )
    pairs = [t for t in tables if "expected_kinship" in t["columns"]][-1]
    assert "observed" not in pairs["columns"]
    # The parents' pair (declared unrelated) and the cross-family one are hidden.
    assert {r["pair"] for r in pairs["rows"]} == {"S1 – S3", "S2 – S3"}
    await user.should_see("not declared related in the pedigree are not shown")


async def test_karyotype_is_neutral_for_readers(user, monkeypatch):
    tables = await open_family_qc(user, monkeypatch, "reader", fid="AU-X-1")
    sex = [t for t in tables if "karyotype" in t["columns"]][-1]
    assert {r["sample"]: r["karyotype"] for r in sex["rows"]} == {
        "S5": "Atypical",
        "S6": "XY",
    }


async def test_family_header_badge_and_member_icons(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("QC: 1 high · 1 medium", retries=10)
    icons = [i for i in user.find(ui.icon).elements if i.props.get("name") == "warning"]
    # S1 and S3: NOT_PARENT is about both; S3's contamination shares its icon.
    assert len(icons) == 2
    # Clicking the badge opens the tab.
    user.find("QC: 1 high · 1 medium").trigger("click")
    await user.should_see("Sex and karyotype", retries=30)


async def test_family_without_qc_has_the_tab_disabled(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1002")
    await user.should_see("Family: 1002", retries=10)
    assert "disable" in _tab(user, "Sample QC").props
    await user.should_not_see("QC ok")


async def test_standalone_family_qc_has_no_cohort_context(user):
    await login(user, "curator")
    await user.open("/family/1001")
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    await user.should_see("open the family from its cohort", retries=10)


async def test_denovo_config_carries_the_qc_warning(user):
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    await user.should_not_see("de novo calls of this family may be unreliable")
    open_tab(user, "wombat-subtabs", "rare_denovo")
    await user.should_see("de novo calls of this family may be unreliable", retries=30)
    await user.should_see("Not the declared parent", retries=10)


async def test_cohort_qc_panel_follows_the_filters(user, monkeypatch):
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-panel-tabs", "Sample QC")
    await user.should_see("Flagged individuals", retries=30)
    await user.should_see("5 of 6 shown individuals checked", retries=30)
    flagged = [t for t in tables if "issues" in t["columns"]][-1]
    # Most severe first, then the most flags: two each for S3 and S5.
    assert [r["sample"] for r in flagged["rows"]][:2] == ["S3", "S5"]
    assert {r["sample"] for r in flagged["rows"]} == {"S1", "S2", "S3", "S5", "S6"}
    await user.should_see("1 related pair(s) across families", retries=10)

    individuals = tables[0]
    with user.client:
        individuals["on_filter"]({"filters": {"FID": "AU"}})
    await user.should_see("2 of 2 shown individuals checked", retries=30)
    flagged = [t for t in tables if "issues" in t["columns"]][-1]
    assert {r["sample"] for r in flagged["rows"]} == {"S5", "S6"}


async def test_cohort_table_qc_column_filters(user, monkeypatch):
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    individuals = tables[0]
    assert "QC" in individuals["columns"]
    by_sample = {r["Sample ID"]: r for r in individuals["rows"]}
    assert by_sample["S4"]["_qc_cats"] == ["Not checked"]
    assert [b["label"] for b in by_sample["S5"]["QC_badges"]] == [
        "Across families",
        "Sex chromosomes",
    ]
    with user.client:
        individuals["on_filter"]({"filters": {"QC": ["Sex chromosomes"]}})
    await user.should_see("1 individuals (of 6 total)", retries=10)
    with user.client:
        individuals["on_filter"]({"filters": {"QC": ["Any issue"]}})
    await user.should_see("5 individuals (of 6 total)", retries=10)


async def test_home_card_shows_the_qc_figure(user):
    await login(user, "reader")
    await user.should_see("QC issues", retries=30)
