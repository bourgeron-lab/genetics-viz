"""The family report buttons, dialog, conclusion and print, on both family pages."""

import asyncio
import sys

from ui_helpers import login


async def _frame(user):
    """The report iframe once the off-thread build is in."""
    await user.should_see(marker="family-report-frame", retries=50)
    return user.find(marker="family-report-frame").elements.pop()


async def test_curator_opens_the_report(user, data_dir):
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    user.find(marker="family-report-view").click()
    await user.should_see(marker="family-report-dialog", retries=10)
    frame = await _frame(user)
    html = frame.props["srcdoc"]
    assert "Genomic analysis report" in html and "Family 1001" in html
    assert "Not the declared parent" in html  # curator wording
    assert frame.props["title"].startswith("DEMO_1001_report_")


async def test_viewer_has_no_report_buttons(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    await user.should_not_see(marker="family-report-view")
    await user.should_not_see(marker="family-report-pdf")


async def test_standalone_family_page_has_the_report(user):
    await login(user, "curator")
    await user.open("/family/AU-X-1")
    await user.should_see("Family: AU-X-1", retries=10)
    user.find(marker="family-report-view").click()
    html = (await _frame(user)).props["srcdoc"]
    assert "Relatedness with other families was not checked" in html
    assert "XXY" in html


async def test_conclusion_is_saved_and_printed(user, data_dir):
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    user.find(marker="family-report-view").click()
    await _frame(user)
    user.find(
        marker="family-report-summary"
    ).elements.pop().value = "Recontact advised."
    user.find(marker="family-report-save").click()
    await user.should_see("Conclusion saved", retries=30)
    rows = (data_dir / "reports" / "family_summaries.tsv").read_text().splitlines()
    assert rows[0] == "FID\tSummary\tUser\tTimestamp"
    assert rows[1].split("\t")[:3] == ["1001", "Recontact advised.", "curator"]
    # The save rebuilds the report: the frame goes, then comes back with it.
    html = ""
    for _ in range(100):
        frames = list(user.client.layout.descendants())
        html = next(
            (
                e.props.get("srcdoc", "")
                for e in frames
                if "family-report-frame" in (e._markers or [])
            ),
            "",
        )
        if "Recontact advised." in html:
            break
        await asyncio.sleep(0.05)
    assert "Recontact advised." in html


async def test_pdf_button_prints_the_report(user, monkeypatch):
    await login(user, "curator")
    module = sys.modules["genetics_viz.components.family_report_dialog"]
    printed = []
    monkeypatch.setattr(module, "print_html", lambda html, title: printed.append(title))
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    user.find(marker="family-report-pdf").click()
    for _ in range(100):
        if printed:
            break
        await asyncio.sleep(0.05)
    assert printed and printed[0].startswith("DEMO_1001_report_")
