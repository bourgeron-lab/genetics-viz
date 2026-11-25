"""Aneuploidy review: the SVs tab's card and the Sample QC karyotype flag."""

import datagen
from nicegui import ui
from ui_helpers import capture_tables, login, open_tab


def _close_dialogs(user):
    for d in user.find(ui.dialog).elements:
        d.close()


async def _open_family(user, monkeypatch, username, fid):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open(f"/cohort/DEMO/family/{fid}")
    await user.should_see(f"Family: {fid}", retries=10)
    return tables


async def test_curator_validates_a_trisomy_from_the_svs_tab(
    user, monkeypatch, data_dir
):
    datagen.add_report_case(data_dir)  # S3: a chr21 gain in WisecondorX
    await _open_family(user, monkeypatch, "curator", "1001")
    open_tab(user, "family-tabs", "SVs")
    await user.should_see(marker="aneuploidy-card", retries=30)
    await user.should_see("1 to review")
    await user.should_see(marker="aneuploidy-row-S3-21")
    user.find(marker="aneuploidy-review-S3-21").click()
    await user.should_see("Aneuploidy — S3, chr21 (+21)", retries=10)
    await user.should_see(marker="aneuploidy-wcx-evidence")
    await user.should_see(marker="aneuploidy-segments")
    karyotype = user.find(marker="aneuploidy-karyotype").elements.pop()
    assert karyotype.value == "47,XY,+21"
    assert user.find(marker="aneuploidy-direction").elements.pop().value == "gain"
    user.find(marker="aneuploidy-origin").elements.pop().value = "de novo"
    user.find(marker="aneuploidy-validation-save").click()
    await user.should_see("Aneuploidy validation saved", retries=30)

    header, row = (
        (data_dir / "validations" / "aneuploidies.tsv").read_text().splitlines()
    )
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert {k: fields[k] for k in ("FID", "Variant", "Sample", "Validation")} == {
        "FID": "1001",
        "Variant": "chr21:ANEUPLOIDY",
        "Sample": "S3",
        "Validation": "present",
    }
    assert (fields["Inheritance"], fields["Direction"], fields["Karyotype"]) == (
        "de novo",
        "gain",
        "47,XY,+21",
    )

    _close_dialogs(user)
    await user.should_not_see("1 to review", retries=30)
    await user.should_see("present")


async def test_curator_validates_a_karyotype_from_the_qc_issues(
    user, monkeypatch, data_dir
):
    tables = await _open_family(user, monkeypatch, "curator", "AU-X-1")
    # S6's sex mismatch is acknowledged in the test data: S5's XXY is left.
    await user.should_see("QC: 1 high", retries=10)
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    await user.should_not_see("Reading the cohort's tables", retries=30)
    issues = [t for t in tables if "severity" in t["columns"]][-1]
    row = next(r for r in issues["rows"] if r["flag"] == "Karyotype XXY")
    with user.client:
        issues["on_row_action"]({"action": "review_flag", "row": row})
    await user.should_see("Aneuploidy — S5, sex chromosomes (XXY)", retries=10)
    await user.should_see(marker="aneuploidy-qc-evidence")
    await user.should_see(marker="aneuploidy-flag-evidence")
    assert user.find(marker="aneuploidy-karyotype").elements.pop().value == "47,XXY"
    user.find(marker="aneuploidy-validation-save").click()
    await user.should_see("Aneuploidy validation saved", retries=30)
    _close_dialogs(user)
    # The validation acknowledges the karyotype flag: the header follows.
    await user.should_see("QC ok", retries=30)


async def test_reader_sees_the_card_without_the_form(user, monkeypatch, data_dir):
    datagen.add_report_case(data_dir)
    await _open_family(user, monkeypatch, "reader", "1001")
    open_tab(user, "family-tabs", "SVs")
    await user.should_see(marker="aneuploidy-row-S3-21", retries=30)
    await user.should_not_see("1 to review")
    user.find(marker="aneuploidy-review-S3-21").click()
    await user.should_see("Aneuploidy — S3, chr21 (+21)", retries=10)
    await user.should_see("Curators and administrators can record validations.")
    await user.should_not_see(marker="aneuploidy-validation-save")
