"""UPD validation and diagnostic from the ROH and Sample QC tabs, and the
family Diagnostics panel's STR and UPD lines."""

import asyncio

import datagen
from nicegui import ui
from ui_helpers import badge_becomes, badge_text, capture_tables, login, open_tab


async def _open_family(user, monkeypatch, username):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    return tables


async def _open_roh_upd(user):
    open_tab(user, "family-tabs", "ROH")
    await user.should_see("UPD candidates", retries=30)
    await user.should_see("S3 · chr2")
    user.find(marker="roh-upd-review-S3").click()
    await user.should_see("Uniparental disomy — S3, chr2", retries=10)


def _close_dialogs(user):
    for d in user.find(ui.dialog).elements:
        d.close()


async def test_curator_validates_a_upd_from_the_roh_tab(user, monkeypatch, data_dir):
    datagen.add_upd_case(data_dir)
    await _open_family(user, monkeypatch, "curator")
    await user.should_see("QC: 2 high · 1 medium", retries=10)
    await _open_roh_upd(user)
    await user.should_see("ROH: an isolated block of long runs of homozygosity")
    # The flag names the mother as the parent without an allele: paternal.
    origin = user.find(marker="upd-origin").elements.pop()
    assert origin.value == "paternal"
    user.find(marker="upd-validation-save").click()
    await user.should_see("UPD validation saved", retries=30)

    rows = (data_dir / "validations" / "upds.tsv").read_text().splitlines()
    assert rows[0].split("\t")[:6] == [
        "FID", "Variant", "Sample", "User", "Inheritance", "Validation",
    ]  # fmt: skip
    assert rows[1].split("\t")[:6] == [
        "1001", "chr2:UPD", "S3", "curator", "paternal", "present",
    ]  # fmt: skip

    # Closing the dialog: the UPD acknowledges the QC flag, the header follows.
    _close_dialogs(user)
    await user.should_see("QC: 1 high · 1 medium", retries=30)
    await user.should_see("present · paternal", retries=10)


async def test_curator_validates_a_upd_from_the_qc_issues(user, monkeypatch, data_dir):
    datagen.add_upd_case(data_dir)
    tables = await _open_family(user, monkeypatch, "curator")
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    await user.should_not_see("Reading the cohort's tables", retries=30)
    issues = [t for t in tables if "severity" in t["columns"]][-1]
    row = next(r for r in issues["rows"] if r["flag"] == "Uniparental disomy candidate")
    with user.client:
        issues["on_row_action"]({"action": "review_flag", "row": row})
    await user.should_see("Uniparental disomy — S3, chr2", retries=10)
    await user.should_see(marker="upd-qc-evidence")
    user.find(marker="upd-validation").elements.pop().value = "absent"
    built = len(tables)
    user.find(marker="upd-validation-save").click()
    await user.should_see("UPD validation saved", retries=30)
    _close_dialogs(user)
    await user.should_see("QC: 1 high · 1 medium", retries=30)
    for _ in range(100):
        if len(tables) > built:
            break
        await asyncio.sleep(0.05)
    issues = [t for t in tables if "severity" in t["columns"]][-1]
    upd = next(r for r in issues["rows"] if r["flag"] == "Uniparental disomy candidate")
    assert upd["review"] == "Explained"
    assert upd["note"] == "UPD validated absent (paternal)"


async def test_reader_sees_no_upd_form(user, monkeypatch, data_dir):
    datagen.add_upd_case(data_dir)
    await _open_family(user, monkeypatch, "reader")
    await _open_roh_upd(user)
    await user.should_see("Curators and administrators can record validations.")
    await user.should_not_see(marker="upd-validation-save")
    # The flag in reader wording, without the pipeline's evidence.
    await user.should_see("Sample QC: Uniparental disomy (S2 – S3)")
    await user.should_not_see("opposite homozygotes")


def _write_diagnostics(data_dir):
    diagnostics = data_dir / "diagnostics"
    diagnostics.mkdir(parents=True, exist_ok=True)
    header = (
        "FID\tVariant\tGene\tImpact\tSample\tUser\tTimestamp\tComment\tIgnore"
        "\tDiagnostic\n"
    )
    (diagnostics / "strs.tsv").write_text(
        header + "1001\tchr4:3074876-3074933:HTT\tHTT\tpathogenic\tS3\tcurator"
        "\t2026-10-12T10:00:00\t\t0\tpathogenic\n"
    )
    (diagnostics / "upds.tsv").write_text(
        header + "1001\tchr2:UPD\t\tROH + Sample QC\tS3\tcurator"
        "\t2026-10-12T10:00:00\t\t0\tuncertain\n"
    )


async def test_diagnostics_panel_reads_str_and_upd(user, monkeypatch, data_dir):
    _write_diagnostics(data_dir)
    await _open_family(user, monkeypatch, "reader")
    await user.should_see("~40 repeats", retries=10)
    await user.should_see(marker="diagnostic-str")
    await user.should_not_see("chr4:3074876-3074933:HTT")
    await user.should_see(marker="diagnostic-upd")
    await user.should_see("chr2")


async def test_diagnostics_panel_prefers_curated_repeats(user, monkeypatch, data_dir):
    _write_diagnostics(data_dir)
    (data_dir / "validations" / "strs.tsv").write_text(
        "FID\tVariant\tSample\tUser\tInheritance\tValidation\tCuratedRepeats"
        "\tComment\tIgnore\tTimestamp\n"
        "1001\tchr4:3074876-3074933:HTT\tS3\tcurator\tde novo\tpresent\t41"
        "\tlong-read\t0\t2026-10-12T10:00:00\n"
    )
    await _open_family(user, monkeypatch, "reader")
    await user.should_see("41 repeats (curated)", retries=10)


async def test_diagnostic_pages_list_upd(user, monkeypatch, data_dir):
    _write_diagnostics(data_dir)
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/diagnostic/all")
    await user.should_see("2 diagnostics", retries=30)
    row = next(r for r in tables[-1]["rows"] if r["Type"] == "UPD")
    assert row["ImpactBadges"][0]["label"] == "ROH + Sample QC"
    await user.should_not_see("Type: SNV, SV, STR, UPD")


# ---------------------------------------------------------------------------
# The UPD tab (the upd step's tables, datagen.add_upd_calls)
# ---------------------------------------------------------------------------


def _tab_labels(user):
    """The family tab bar's labels, in order (its tooltips sit between them)."""
    bar = user.find(marker="family-tabs").elements.pop()
    return [
        c.props.get("label") for c in bar.default_slot.children if isinstance(c, ui.tab)
    ]


def _candidates(tables):
    return [t for t in tables if "sources" in t["columns"]][-1]


async def _open_upd_tab(user, monkeypatch, username="curator"):
    tables = await _open_family(user, monkeypatch, username)
    open_tab(user, "family-tabs", "UPD")
    await user.should_see(marker="upd-candidates", retries=30)
    return tables


async def test_upd_tab_sits_between_roh_and_str(user, monkeypatch, data_dir):
    datagen.add_upd_calls(data_dir)
    await _open_family(user, monkeypatch, "curator")
    labels = _tab_labels(user)
    assert labels[labels.index("ROH") + 1] == "UPD"
    assert labels[labels.index("UPD") + 1] == "STR"
    # S2's medium and S3's high event wait for a curator.
    assert await badge_text(user, "upd-tab-count") == "2"


async def test_readers_get_no_upd_count(user, monkeypatch, data_dir):
    datagen.add_upd_calls(data_dir)
    await _open_family(user, monkeypatch, "reader")
    await user.should_see(marker="qc-tab-count", retries=50)
    await user.should_not_see(marker="upd-tab-count")


async def test_standalone_family_has_the_upd_tab(user, data_dir):
    datagen.add_upd_calls(data_dir)
    await login(user, "reader")
    await user.open("/family/1001")
    await user.should_see("1001", retries=10)
    labels = _tab_labels(user)
    assert labels[labels.index("ROH") + 1] == "UPD"


async def test_upd_tab_without_any_upd_data_is_disabled(user, monkeypatch, data_dir):
    await _open_family(user, monkeypatch, "reader")
    bar = user.find(marker="family-tabs").elements.pop()
    upd = next(c for c in bar.default_slot.children if c.props.get("label") == "UPD")
    assert "disable" in upd.props


async def test_curator_reviews_an_event_in_the_upd_tab(user, monkeypatch, data_dir):
    datagen.add_upd_calls(data_dir)
    tables = await _open_upd_tab(user, monkeypatch)
    await user.should_see("upd step · upd_v4")
    candidates = _candidates(tables)
    assert [(r["sample"], r["chrom"]) for r in candidates["rows"]] == [
        ("S2", "chr15"),
        ("S3", "chr2"),
    ]
    s3 = candidates["rows"][1]
    assert s3["sources_badges"][0]["label"] == "UPD · High"
    assert s3["type"] == "isodisomy · segmental · terminal"
    assert s3["validation"] == "not reviewed"
    await user.should_see(marker="upd-genome-svg", retries=30)

    with user.client:
        candidates["on_row_action"]({"action": "review_upd", "row": s3})
    await user.should_see("Uniparental disomy — S3, chr2", retries=10)
    await user.should_see(marker="upd-step-evidence")
    await user.should_see(marker="upd-chromosome-chart", retries=30)
    chart = user.find(marker="upd-chromosome-chart").elements.pop()
    names = [s["name"] for s in chart.options["series"]]
    assert "Lacks the mother's allele" in names and "Het rate" in names
    # The form starts from the event.
    assert user.find(marker="upd-origin").elements.pop().value == "paternal"
    assert user.find(marker="upd-composition").elements.pop().value == "iso"
    assert user.find(marker="upd-start").elements.pop().value == 0.1
    user.find(marker="upd-validation-save").click()
    await user.should_see("UPD validation saved", retries=30)

    row = upd_rows(data_dir)[-1]
    assert (row["Variant"], row["Composition"], row["Extent"], row["UpdLabel"]) == (
        "chr2:UPD", "iso", "segmental", "upd_v4",
    )  # fmt: skip
    assert (row["Start"], row["End"]) == ("100000", "39700000")

    built = len(tables)
    _close_dialogs(user)
    for _ in range(100):
        if len(tables) > built:
            break
        await asyncio.sleep(0.05)
    s3 = next(r for r in _candidates(tables)["rows"] if r["sample"] == "S3")
    assert s3["validation"] == "present"
    assert await badge_becomes(user, "upd-tab-count", "1")


def upd_rows(data_dir):
    import csv

    with open(data_dir / "validations" / "upds.tsv") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


async def test_upd_tab_shows_low_events_and_older_versions(user, monkeypatch, data_dir):
    datagen.add_upd_calls(data_dir)
    tables = await _open_upd_tab(user, monkeypatch, "reader")
    user.find(marker="upd-confidence-low").elements.pop().value = True
    for _ in range(100):
        rows = _candidates(tables)["rows"]
        if len(rows) == 3:
            break
        await asyncio.sleep(0.05)
    assert ("S3", "chr16") in {(r["sample"], r["chrom"]) for r in rows}
    chr16 = next(r for r in rows if r["chrom"] == "chr16")
    assert chr16["validation"] == "not reviewed"

    user.find(marker="upd-label").elements.pop().value = "upd_v3"
    await user.should_see("upd step · upd_v3", retries=10)
    await user.should_see("wrote no per-window genotype counts", retries=10)


async def test_reader_sees_neutral_members_and_no_form(user, monkeypatch, data_dir):
    datagen.add_upd_calls(data_dir)
    tables = await _open_upd_tab(user, monkeypatch, "reader")
    members = user.find(marker="upd-members").elements.pop()
    table = next(e for e in members.descendants() if isinstance(e, ui.table))
    assert "parents" not in {c["name"] for c in table.columns}
    candidates = _candidates(tables)
    with user.client:
        candidates["on_row_action"](
            {"action": "review_upd", "row": candidates["rows"][1]}
        )
    await user.should_see("Uniparental disomy — S3, chr2", retries=10)
    await user.should_see("Curators and administrators can record validations.")
    await user.should_not_see(marker="upd-validation-save")


async def test_a_link_opens_the_review(user, data_dir):
    datagen.add_upd_calls(data_dir)
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001?tab=upd&upd=S3:2")
    await user.should_see("Uniparental disomy — S3, chr2", retries=30)


async def test_an_unsettled_qc_flag_is_reviewed_in_the_upd_tab(
    user, monkeypatch, data_dir
):
    datagen.add_upd_case(data_dir)
    datagen.add_upd_calls(data_dir)
    tables = await _open_family(user, monkeypatch, "curator")
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    await user.should_not_see("Reading the cohort's tables", retries=30)
    issues = [t for t in tables if "severity" in t["columns"]][-1]
    row = next(r for r in issues["rows"] if r["flag"] == "Uniparental disomy candidate")
    with user.client:
        issues["on_row_action"]({"action": "review_flag", "row": row})
    # The UPD tab opens, with the review and its chromosome view.
    await user.should_see(marker="upd-candidates", retries=30)
    await user.should_see("Uniparental disomy — S3, chr2", retries=10)
    await user.should_see(marker="upd-chromosome-chart", retries=30)
    bar = user.find(marker="family-tabs").elements.pop()
    assert bar.value.props.get("name") == "UPD"
    s3 = next(r for r in _candidates(tables)["rows"] if r["sample"] == "S3")
    assert [b["label"] for b in s3["sources_badges"]] == [
        "UPD · High",
        "ROH",
        "Sample QC",
    ]


async def test_downloads_follow_the_role(user, data_dir):
    import csv
    import io

    datagen.add_upd_calls(data_dir)
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    open_tab(user, "family-tabs", "UPD")
    await user.should_see(marker="upd-candidates", retries=30)
    user.find(marker="download-upd-members").click()
    response = await user.download.next(timeout=5)
    reader = csv.DictReader(io.StringIO(response.content.decode()), delimiter="\t")
    assert "father_used" not in reader.fieldnames
    assert [r["IID"] for r in reader] == ["S1", "S2", "S3"]
    user.find(marker="download-upd-windows").click()
    response = await user.download.next(timeout=5)
    assert response.content.decode().startswith("FID\tIID\tmode\tchrom")


# ---------------------------------------------------------------------------
# The cohort page
# ---------------------------------------------------------------------------


async def test_cohort_upd_panel_and_column(user, monkeypatch, data_dir):
    datagen.add_upd_calls(data_dir)
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    individuals = next(t for t in tables if "UPD" in t["columns"])
    by_sample = {r["Sample ID"]: r for r in individuals["rows"]}
    assert [b["label"] for b in by_sample["S3"]["UPD_badges"]] == ["chr2 · High"]
    assert by_sample["S4"]["_upd_cats"] == ["Not analysed"]

    open_tab(user, "cohort-tabs", "UPD")
    await user.should_see(marker="cohort-upd-events", retries=30)
    events = [t for t in tables if "confidence" in t["columns"]][-1]
    assert {(r["sample"], r["chrom"]) for r in events["rows"]} == {
        ("S3", "chr2"),
        ("S2", "chr15"),
    }
    s3 = next(r for r in events["rows"] if r["sample"] == "S3")
    assert s3["_href"] == "/cohort/DEMO/family/1001?tab=upd&upd=S3%3A2"
    chart = user.find(marker="cohort-upd-chart").elements.pop()
    assert chart.options["series"][0]["name"] == "High"

    # The panel follows the table: High keeps S3 only.
    with user.client:
        individuals["on_filter"]({"filters": {"UPD": ["High"]}})
    await user.should_see("1 of 6 individuals", retries=10)
    for _ in range(100):
        events = [t for t in tables if "confidence" in t["columns"]][-1]
        if len(events["rows"]) == 1:
            break
        await asyncio.sleep(0.05)
    assert [r["sample"] for r in events["rows"]] == ["S3"]


async def test_cohort_qc_flag_links_to_the_upd_tab(user, monkeypatch, data_dir):
    datagen.add_upd_case(data_dir)
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-tabs", "Sample QC")
    for _ in range(200):
        flagged = [t for t in tables if "issues" in t["columns"]]
        if flagged:
            break
        await asyncio.sleep(0.05)
    s3 = next(r for r in flagged[-1]["rows"] if r["sample"] == "S3")
    assert s3["_href"] == "/cohort/DEMO/family/1001?tab=upd&upd=S3%3A2"
    s1 = next(r for r in flagged[-1]["rows"] if r["sample"] != "S3")
    assert s1["_href"] is None


async def test_cell_line_caveat_on_the_upd_tab_and_dialog(user, monkeypatch, data_dir):
    # S3 has a chr2 UPD and its DNA is from a cell line (genome QC).
    datagen.add_upd_calls(data_dir)
    datagen.add_genomeqc_tables(data_dir)
    tables = await _open_upd_tab(user, monkeypatch)
    await user.should_see(marker="cell-line-note-upd")
    await user.should_see("S3: DNA from a cell line.")
    s3 = next(r for r in _candidates(tables)["rows"] if r["sample"] == "S3")
    with user.client:
        _candidates(tables)["on_row_action"]({"action": "review_upd", "row": s3})
    await user.should_see("Uniparental disomy — S3, chr2", retries=10)
    notes = user.find(marker="cell-line-note-upd").elements
    assert len(notes) == 2  # the tab's and the dialog's
