"""UPD validation and diagnostic from the ROH and Sample QC tabs, and the
family Diagnostics panel's STR and UPD lines."""

import asyncio

import datagen
from nicegui import ui
from ui_helpers import capture_tables, login, open_tab


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
