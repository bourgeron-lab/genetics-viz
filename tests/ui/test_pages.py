"""Pages rendered in NiceGUI's user simulation, logged in as each role."""

import pytest
from nicegui import ui
from ui_helpers import capture_tables, login, open_tab, roh_threshold_labels


async def test_health(user):
    response = await user.http_client.get("/health")
    assert response.status_code == 200 and response.json()["status"] == "ok"


async def test_anonymous_visitor_is_sent_to_login(user):
    await user.open("/")
    await user.should_see("Sign in to continue")


async def test_wrong_password_is_refused(user):
    await user.open("/login")
    user.find("Username").type("reader")
    user.find("Password").type("not-the-password")
    user.find(kind=ui.button, content="Sign in").click()
    await user.should_see("Invalid username or password")


async def test_home_lists_complete_and_incomplete_cohorts(user):
    await login(user, "reader")
    await user.should_see("DEMO")
    await user.should_see("NOHEAD")
    await user.should_see("EMPTY")
    await user.should_see("BROKEN")


async def test_unknown_cohort_is_not_found(user):
    await login(user, "reader")
    await user.open("/cohort/NOPE")
    await user.should_see("not found", retries=10)


async def test_cohort_page_and_its_tabs(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    for label, expected in [
        ("Ancestry", "individuals"),
        ("Polygenic Scores", "Polygenic Scores"),
        ("ROH", "Runs of homozygosity"),
    ]:
        open_tab(user, "cohort-panel-tabs", label)
        await user.should_see(expected, retries=30)
    plots = [e for e in user.find(ui.plotly).elements]
    assert plots, "the ROH and PGS tabs draw Plotly violins"


async def test_family_page_and_its_tabs(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    open_tab(user, "family-tabs", "Ancestry")
    await user.should_see("bundle 1.1.0 / dp10gq20", retries=30)
    open_tab(user, "family-tabs", "Polygenic Scores")
    await user.should_see("2 polygenic scores", retries=30)
    open_tab(user, "family-tabs", "ROH")
    await user.should_see("Runs of homozygosity", retries=30)
    await user.should_see("Compared with the reference", retries=30)
    charts = user.find(ui.echart).elements
    assert any("reference" in str(c.options.get("title", {})).lower() for c in charts)


async def open_family_roh(user, monkeypatch, username):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    open_tab(user, "family-tabs", "ROH")
    await user.should_see("Compared with the reference", retries=30)
    return tables


async def test_roh_relationship_wording_is_hidden_from_readers(user, monkeypatch):
    tables = await open_family_roh(user, monkeypatch, "reader")
    summary = next(t for t in tables if "froh" in t["columns"])
    assert "parents_kinship" not in summary["columns"]
    assert all(
        "degree" not in b["label"] for r in summary["rows"] for b in r["grade_badges"]
    )
    kinship = next(t for t in tables if "kinship" in t["columns"])
    assert {r["declared"] for r in kinship["rows"]} == {"Parent–offspring"}
    assert not any(
        "cousin" in label or "degree" in label for label in roh_threshold_labels(user)
    )


async def test_roh_relationship_wording_is_shown_to_curators(user, monkeypatch):
    tables = await open_family_roh(user, monkeypatch, "curator")
    summary = next(t for t in tables if "froh" in t["columns"])
    assert "parents_kinship" in summary["columns"]
    child = next(r for r in summary["rows"] if r["IID"] == "S3")
    assert child["grade_badges"][0]["label"] == "High · second-degree-like"
    kinship = next(t for t in tables if "kinship" in t["columns"])
    assert "Parents of a child" in {r["declared"] for r in kinship["rows"]}
    assert any("first-cousin-like" in label for label in roh_threshold_labels(user))


async def test_workflow_tab_shows_the_run_record_and_parameters(user):
    import sys
    from datetime import datetime

    docs = sys.modules["genetics_viz.utils.ghfc_ngs_docs"]
    # A fresh cached document: the tab reads it instead of fetching GitHub.
    docs._cached = docs.FetchedDoc(
        markdown="# Stub reference", error=None, fetched_at=datetime.now()
    )
    await login(user, "reader")
    await user.open("/cohort/DEMO")
    open_tab(user, "cohort-page-tabs", "Workflow")
    await user.should_see("Runs of homozygosity", retries=30)
    await user.should_see("Ancestry & PGS", retries=30)
    await user.should_see("not measured", retries=30)


def _tab(user, name):
    return next(t for t in user.find(ui.tab).elements if t.props.get("name") == name)


async def test_family_str_tab_pivots_members_and_adds_cohort_context(user, monkeypatch):
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    open_tab(user, "family-tabs", "STR")
    await user.should_see("ExpansionHunter", retries=30)
    open_tab(user, "str-tabs", "EH")
    await user.should_see("compared with 3 founders", retries=30)

    table = tables[-1]
    assert {"m0", "m1", "m2", "founders_max", "cohort_flagged"} <= set(table["columns"])
    first = table["rows"][0]
    assert first["locus"] == "HTT"
    assert [b["label"] for b in first["m2_badges"]] == ["18", "40·IR"]
    assert _tab(user, "EH").props.get("label") == "EH (2 flagged)"


async def test_family_without_str_has_the_tab_disabled(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1002")
    await user.should_see("Family: 1002", retries=10)
    assert "disable" in _tab(user, "STR").props


async def test_standalone_family_str_tab_has_no_cohort_context(user):
    await login(user, "reader")
    await user.open("/family/AU-X-1")
    open_tab(user, "family-tabs", "STR")
    await user.should_see("ExpansionHunter", retries=30)
    await user.should_see("No cohort context on this page.", retries=30)


@pytest.mark.xfail(
    strict=True,
    reason="load_family_pedigree reads the first row of the pipeline's headerless"
    " family pedigree as a header, so the standalone page finds no members",
)
async def test_standalone_family_reads_a_headerless_family_pedigree(user):
    await login(user, "reader")
    await user.open("/family/AU-X-1")
    await user.should_see("Family: AU-X-1", retries=10)
    await user.should_not_see("No pedigree file found for this family.")
