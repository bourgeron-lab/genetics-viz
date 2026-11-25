"""The diagnostic dialog's five tiers and ACMG criteria checklist (#44)."""

from ui_helpers import capture_tables, login

from genetics_viz.components.diagnostic_loader import DIAGNOSTIC_HEADER


async def _open_shank3_dialog(user, monkeypatch, username):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    shank3 = next(r for r in tables[0]["rows"] if r["VEP_SYMBOL"] == "SHANK3")
    sent = {k: v for k, v in shank3.items() if not k.startswith("_")}
    with user.client:
        tables[0]["on_row_action"]({"action": "open_diagnostic", "row": sent})
    await user.should_see("Diagnostic Review", retries=30)


async def test_criteria_set_the_class_and_are_saved(user, monkeypatch, data_dir):
    await _open_shank3_dialog(user, monkeypatch, "curator")
    verdict = user.find(marker="diagnostic-verdict").elements.pop()
    assert verdict.value == "uncertain"
    assert list(verdict.options) == [
        "pathogenic", "likely pathogenic", "uncertain", "likely benign", "benign"
    ]  # fmt: skip
    # De novo (confirmed parentage) and absent from gnomAD (proposed at
    # Supporting): 4 + 1 = 5 points, uncertain.
    user.find(marker="acmg-PS2").click()
    user.find(marker="acmg-PM2").click()
    await user.should_see("5 points")
    assert verdict.value == "uncertain"
    # PP3 at Moderate: 7 points, likely pathogenic.
    user.find(marker="acmg-PP3").click()
    await user.should_see(marker="acmg-PP3-strength")
    user.find(marker="acmg-PP3-strength").click()
    user.find(marker="acmg-PP3-Moderate").click()
    await user.should_see("7 points")
    assert verdict.value == "likely pathogenic"
    await user.should_not_see(marker="acmg-override")
    # The curator records pathogenic: an override, said before the save.
    verdict.value = "pathogenic"
    await user.should_see(marker="acmg-override")
    user.find(marker="diagnostic-save").click()
    await user.should_see("Diagnostic saved: Pathogenic", retries=30)

    header, line = (data_dir / "diagnostics" / "snvs.tsv").read_text().splitlines()
    assert header.split("\t") == DIAGNOSTIC_HEADER
    saved = dict(zip(DIAGNOSTIC_HEADER, line.split("\t")))
    assert (saved["Diagnostic"], saved["Criteria"], saved["Points"]) == (
        "pathogenic",
        "PS2;PM2_Supporting;PP3_Moderate",
        "7",
    )
    assert saved["Ruleset"] == "richards2015+tavtigian2020"


async def test_saving_without_criteria_writes_none(user, monkeypatch, data_dir):
    await _open_shank3_dialog(user, monkeypatch, "curator")
    user.find(marker="diagnostic-verdict").elements.pop().value = "likely benign"
    user.find(marker="diagnostic-save").click()
    await user.should_see("Diagnostic saved: Likely benign", retries=30)
    _, line = (data_dir / "diagnostics" / "snvs.tsv").read_text().splitlines()
    saved = dict(zip(DIAGNOSTIC_HEADER, line.split("\t")))
    assert (saved["Diagnostic"], saved["Criteria"], saved["Points"]) == (
        "likely benign",
        "",
        "",
    )


async def test_readers_get_no_checklist(user, monkeypatch):
    await _open_shank3_dialog(user, monkeypatch, "reader")
    await user.should_see("Read-only access")
    await user.should_not_see(marker="acmg-criteria")
