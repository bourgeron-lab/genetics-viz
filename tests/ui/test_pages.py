"""Pages rendered in NiceGUI's user simulation, logged in as each role."""

import asyncio

import datagen
from nicegui import ui
from nicegui.testing.user_interaction import UserInteraction
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


async def test_wombat_sub_tab_builds_its_table_when_opened(user, monkeypatch):
    """Each Wombat config's table is built when its sub-tab is first shown.

    Quasar mounts only the shown panel, so a table built at page load for a
    hidden sub-tab had no container in the browser until the tab was opened.
    """
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    assert [len(t["rows"]) for t in tables] == [6]
    await user.should_not_see("Wombat Configuration: strict")

    open_tab(user, "wombat-subtabs", "strict")
    await user.should_see("Wombat Configuration: strict", retries=30)
    assert [len(t["rows"]) for t in tables] == [6, 2]

    # Back and forth builds nothing again.
    open_tab(user, "wombat-subtabs", "default")
    open_tab(user, "wombat-subtabs", "strict")
    await asyncio.sleep(0.3)
    assert len(tables) == 2


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
    genes = next(t for t in tables if "gene_id" in t["columns"])
    [badge] = genes["rows"][0]["GeneBadges"]
    assert badge["label"] == "GENEA"
    assert badge["href"] == "http://genetrek.test/?search=GENEA"


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
    [locus_badge] = first["LocusBadges"]
    assert locus_badge["label"] == "HTT"
    assert locus_badge["href"] == "http://genetrek.test/?search=HTT"
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


async def test_standalone_family_reads_a_headerless_family_pedigree(user):
    await login(user, "reader")
    await user.open("/family/AU-X-1")
    await user.should_see("Family: AU-X-1", retries=10)
    await user.should_not_see("No pedigree file found for this family.")
    await user.should_see("2 members")
    await user.should_see("1 founders")


async def test_cohort_str_panel_follows_the_filters(user, monkeypatch):
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-panel-tabs", "STR")
    await user.should_see("flagged genotypes", retries=30)

    flagged = next(t for t in tables if "alleles" in t["columns"])
    assert {(r["sample"], r["locus"]) for r in flagged["rows"]} == {
        ("S3", "HTT"),
        ("S2", "FXN"),
    }
    assert all(
        r["LocusBadges"][0]["href"].endswith("?search=" + r["locus"].split("_")[0])
        for r in flagged["rows"]
    )

    individuals = tables[0]
    with user.client:
        individuals["on_filter"]({"filters": {"FID": "AU"}})
    await user.should_see("No flagged genotype among the individuals shown", retries=30)
    loci = [t for t in tables if "n_pathogenic" in t["columns"]][-1]
    assert sum(r["n_flagged"] for r in loci["rows"]) == 0


async def open_str_locus(user, monkeypatch, username, locus="HTT"):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    open_tab(user, "family-tabs", "STR")
    await user.should_see("compared with 3 founders", retries=30)
    with user.client:
        tables[-1]["on_row_action"]({"action": "view_locus", "row": {"locus": locus}})
    await user.should_see("REViewer", retries=30)
    return tables


async def test_curator_saves_an_str_validation(user, monkeypatch, data_dir):
    tables = await open_str_locus(user, monkeypatch, "curator")
    await user.should_see("No validation recorded at this locus.")
    built = len(tables)
    user.find(marker="str-validation-save").click()
    await user.should_see("STR validation saved", retries=30)

    lines = (data_dir / "validations" / "strs.tsv").read_text().splitlines()
    assert lines[1].split("\t")[:6] == [
        "1001", "chr4:3074876-3074933:HTT", "S3", "curator", "unknown", "present"
    ]  # fmt: skip
    # The table refresh is fire-and-forget: wait for the rebuilt table.
    for _ in range(100):
        if len(tables) > built:
            break
        await asyncio.sleep(0.05)
    htt = next(r for r in tables[-1]["rows"] if r["locus"] == "HTT")
    assert htt["m2_badges"][-1]["label"] == "✓"


async def test_reader_sees_no_str_validation_form(user, monkeypatch):
    await open_str_locus(user, monkeypatch, "reader")
    await user.should_see("Curators and administrators can record validations.")
    await user.should_not_see(marker="str-validation-save")


def _write_str_validation(data_dir):
    path = data_dir / "validations" / "strs.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "FID\tVariant\tSample\tUser\tInheritance\tValidation\tCuratedRepeats"
        "\tComment\tIgnore\tTimestamp\n"
        "1001\tchr4:3074876-3074933:HTT\tS3\tcurator\tde novo\tpresent\t41"
        "\tlong-read\t0\t2026-10-12T10:00:00\n"
    )


async def test_validation_pages_list_str(user, monkeypatch, data_dir):
    _write_str_validation(data_dir)
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/validation/all")
    await user.should_see("1 validations", retries=30)
    row = tables[-1]["rows"][0]
    assert (row["Type"], row["Sample"], row["Validation"]) == ("STR", "S3", "present")

    await user.open("/validation/statistics")
    await user.should_see("Validation Statistics", retries=30)
    assert any(c.text == "STR" for c in user.find(ui.checkbox).elements)


async def test_reviewer_images_and_why_one_is_missing(user, monkeypatch):
    await open_str_locus(user, monkeypatch, "reader")
    await user.should_see("S3 — REViewer")
    await user.should_see("S1: not flagged, not drawn")


async def test_flagged_genotype_without_reviewer_image(user, monkeypatch):
    await open_str_locus(user, monkeypatch, "reader", locus="FXN")
    await user.should_see("S2: no REViewer image")


async def test_curator_saves_an_str_diagnostic(user, monkeypatch, data_dir):
    tables = await open_str_locus(user, monkeypatch, "curator")
    built = len(tables)
    user.find(marker="str-diagnostic-S3").click()
    await user.should_see("Diagnostic Review", retries=30)
    user.find(marker="diagnostic-save").click()
    await user.should_see("Diagnostic saved", retries=30)

    lines = (data_dir / "diagnostics" / "strs.tsv").read_text().splitlines()
    fields = lines[1].split("\t")
    assert fields[:5] == [
        "1001", "chr4:3074876-3074933:HTT", "HTT", "pathogenic", "S3"
    ]  # fmt: skip
    assert fields[-1] == "uncertain"
    # The section refreshes in place; the family table waits for the close.
    await user.should_see("uncertain", retries=30)
    assert len(tables) == built
    with user.client:
        for dialog in user.find(ui.dialog).elements:
            dialog.close()
    for _ in range(100):
        if len(tables) > built:
            break
        await asyncio.sleep(0.05)
    htt = next(r for r in tables[-1]["rows"] if r["locus"] == "HTT")
    assert htt["m2_badges"][-1]["label"] == "Dx ?"


async def test_reader_cannot_save_an_str_diagnostic(user, monkeypatch):
    await open_str_locus(user, monkeypatch, "reader")
    user.find(marker="str-diagnostic-S3").click()
    await user.should_see("Read-only access", retries=30)
    await user.should_not_see(marker="diagnostic-save")


def _write_str_diagnostic(data_dir):
    path = data_dir / "diagnostics" / "strs.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "FID\tVariant\tGene\tImpact\tSample\tUser\tTimestamp\tComment\tIgnore"
        "\tDiagnostic\n"
        "1001\tchr4:3074876-3074933:HTT\tHTT\tpathogenic\tS3\tcurator"
        "\t2026-10-12T10:00:00\t\t0\tpathogenic\n"
    )


async def test_diagnostic_pages_list_str(user, monkeypatch, data_dir):
    _write_str_diagnostic(data_dir)
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/diagnostic/all")
    await user.should_see("1 diagnostics", retries=30)
    row = tables[-1]["rows"][0]
    assert (row["Type"], row["Sample"], row["Diagnostic"]) == (
        "STR",
        "S3",
        "pathogenic",
    )
    assert row["ImpactBadges"][0]["label"] == "Pathogenic"

    await user.open("/diagnostic/statistics")
    await user.should_see("Diagnostic Statistics", retries=30)
    assert any(c.text == "STR" for c in user.find(ui.checkbox).elements)


async def test_sv_diagnostic_saved_from_search_is_read_back(
    user, monkeypatch, data_dir
):
    """The search page saves an SV diagnostic under the key its table reads."""
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/search/DEMO")
    await user.should_see("Search: DEMO", retries=30)
    user.find(kind=ui.button, content="+ WisecondorX").click()
    await user.should_see("Exonic only", retries=30)
    user.find("Locus (optional)").type("chr2:165000000-165500000")
    user.find(marker="search-button").click()
    for _ in range(100):
        if tables:
            break
        await asyncio.sleep(0.05)
    [row] = tables[-1]["rows"]
    assert (row["sample"], row["Diagnostic"]) == ("S3", "")
    # The SV's region opens GeneTrek on its genes, in the fragment.
    assert row["_genetrek_url"] == (
        f"http://genetrek.test/#search={row['Variant']}&genome=GRCh38"
    )

    built = len(tables)
    with user.client:
        tables[-1]["on_row_action"]({"action": "open_diagnostic", "row": row})
    await user.should_see("Diagnostic Review", retries=30)
    user.find(marker="diagnostic-save").click()

    lines = (data_dir / "diagnostics" / "svs.tsv").read_text().splitlines()
    assert lines[1].split("\t")[:2] == ["1001", "chr2:165114001-165316000:del"]
    for _ in range(100):
        if len(tables) > built:
            break
        await asyncio.sleep(0.05)
    [row] = tables[-1]["rows"]
    assert row["Diagnostic"] == "uncertain"


def _sample_charts(user, marker):
    return [c.options for c in user.find(marker=marker).elements]


async def test_wombat_stats_count_variants_per_sample(user):
    """The Stats dialog counts each member's variants and places them in the cohort."""
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    user.find(kind=ui.button, content="Stats").click()
    await user.should_see("Variants per Sample (current view)", retries=10)

    (bars,) = _sample_charts(user, "sample-count-bars")
    assert bars["xAxis"]["data"] == ["S1", "S2", "S3"]
    totals = [sum(vals) for vals in zip(*(s["data"] for s in bars["series"]))]
    assert totals == [0, 0, 6]

    await user.should_see("Compared with the cohort (unfiltered)")
    await user.should_see("DEMO.rare.vep1.default.results.tsv", retries=10)
    await user.should_see(marker="sample-count-histogram", retries=30)
    histograms = _sample_charts(user, "sample-count-histogram")
    titles = sorted(h["title"]["text"] for h in histograms)
    assert titles == ["Cohort founders (n=4)", "Cohort non-founders (n=2)"]


async def test_sv_stats_count_calls_per_sample(user, data_dir):
    """The WisecondorX Stats dialog does the same against the cohort BED."""
    family_wcx = (
        datagen.family_dir(data_dir, "1001")
        / "svs"
        / "wisecondorx"
        / "1001_aberrations.annotated.bed"
    )
    datagen.write_tsv(family_wcx, datagen.WCX_COLUMNS, datagen.WCX_ROWS[:1])

    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    open_tab(user, "family-tabs", "SVs")
    await user.should_see("WisecondorX Aberrations", retries=30)
    # The Wombat tab's Stats button is built first; the SV one is the newest
    stats_buttons = user.find(kind=ui.button, content="Stats").elements
    UserInteraction(user, {max(stats_buttons, key=lambda b: b.id)}, "Stats").click()
    await user.should_see("SVs per Sample (current view)", retries=10)

    bars = [b for b in _sample_charts(user, "sample-count-bars") if b["series"]]
    assert any(
        b["xAxis"]["data"] == ["S1", "S2", "S3"]
        and [sum(v) for v in zip(*(s["data"] for s in b["series"]))] == [0, 0, 1]
        for b in bars
    )
    await user.should_see("DEMO_aberrations.bed", retries=10)


async def test_wombat_stats_read_a_headerless_cohort_file(user, data_dir):
    """The cohort "strict" table has no header; the family file lends its own."""
    datagen.build_headerless_cohort_wombat(data_dir)
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    open_tab(user, "wombat-subtabs", "strict")
    await user.should_see("Wombat Configuration: strict", retries=30)
    # The default sub-tab's Stats button is built first; strict's is the newest
    stats_buttons = user.find(kind=ui.button, content="Stats").elements
    UserInteraction(user, {max(stats_buttons, key=lambda b: b.id)}, "Stats").click()
    await user.should_see("DEMO.rare.vep1.annotated.strict.results.tsv", retries=10)
    await user.should_see(marker="sample-count-histogram", retries=30)
    await user.should_not_see("no comparison available")
    non_founders = next(
        h
        for h in _sample_charts(user, "sample-count-histogram")
        if "non-founders" in h["title"]["text"]
    )
    (line,) = non_founders["series"][0]["markLine"]["data"]
    assert line["label"]["formatter"] == "S3: 2"
