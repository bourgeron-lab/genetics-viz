"""The family status: chip and dialog on the family pages, Families tab, home card.

DEMO holds 1001 (trio, S3 affected), 1002 (S4, phenotype unknown) and AU-X-1
(S5 affected). ``datagen.add_family_statuses`` makes 1001 solved and 1002 not
analysable; AU-X-1 stays not reviewed.
"""

import asyncio
from datetime import datetime
from pathlib import Path

import datagen
from ui_helpers import capture_tables, login

from genetics_viz.utils import family_status as fs
from genetics_viz.utils import pedigree_edit as pe

DEMO = Path("cohorts") / "DEMO" / "DEMO.pedigree.tsv"


async def wait_for(condition, tries: int = 200) -> None:
    for _ in range(tries):
        if condition():
            return
        await asyncio.sleep(0.05)
    assert condition()


def _text(user, marker: str) -> str:
    found = user.find(marker=marker).elements
    return found.pop().text if found else ""


async def test_a_curator_sets_the_status(user, data_dir):
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/AU-X-1")
    await user.should_see("Family: AU-X-1", retries=20)
    await user.should_see(marker="family-status", retries=20)
    assert _text(user, "family-status") == "Status: Not reviewed"
    user.find(marker="family-status").click()
    await user.should_see(marker="family-status-dialog", retries=20)
    user.find(marker="family-status-select").elements.pop().set_value("unsolved")
    user.find(marker="family-status-comment").type("nothing in the panel genes")
    user.find(marker="family-status-save").click()
    await wait_for(lambda: _text(user, "family-status") == "Status: Unsolved")
    record = fs.family_statuses(data_dir)["AU-X-1"]
    assert (record.status, record.user) == ("unsolved", "curator")
    assert record.comment == "nothing in the panel genes"
    # The comment is in the curator's tooltip.
    await user.should_see("nothing in the panel genes")


async def test_a_reader_sees_the_status_and_cannot_change_it(user, data_dir):
    datagen.add_family_statuses(data_dir)
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1002")
    await user.should_see(marker="family-status", retries=20)
    assert _text(user, "family-status") == "Status: Not analysable"
    # A reader's tooltip has no comment.
    await user.should_not_see("tube mix-up at the centre")
    user.find(marker="family-status").click()
    await asyncio.sleep(0.3)
    await user.should_not_see(marker="family-status-dialog")


async def test_a_pedigree_edit_after_the_status_reopens_the_family(
    user, data_dir, monkeypatch
):
    datagen.add_family_statuses(data_dir)
    monkeypatch.setattr(pe, "can_write", lambda: True)
    source = pe.load_family(data_dir, "1001", [data_dir / DEMO])
    draft = source.draft()
    draft.set_phenotype("S1", "2")
    pe.save_pedigree_edit(
        source, draft, user="alice", comment="father affected", now=datetime.now()
    )
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see(marker="family-status-reason", retries=60)
    assert _text(user, "family-status") == "Status: Solved"
    assert _text(user, "family-status-reason") == "Reopened"
    # A phenotype edit reruns nothing: ready to review again.
    await user.should_see("ready to be reviewed again")
    # Solved with no pathogenic diagnostic standing.
    await user.should_see("no pathogenic diagnostic stands")


async def test_the_families_tab_status_column_and_filter(user, data_dir, monkeypatch):
    datagen.add_family_statuses(data_dir)
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=20)
    families = next(t for t in tables if "Structure" in t["columns"])
    assert families["columns"][0] == "Status"
    table = families["table"]
    await wait_for(lambda: any(r.get("Status") for r in table.rows))
    rows = {r["FID"]: r for r in table.rows}
    assert rows["1001"]["Status"] == "Solved"
    assert rows["1002"]["Status"] == "Not analysable"
    assert rows["AU-X-1"]["Status"] == "Not reviewed"
    # Readers get the status, not the comment.
    assert "SHANK3" not in rows["1001"]["Status_badges"][0]["tooltip"]
    await user.should_see(marker="families-status-summary-yield", retries=20)
    assert (
        _text(user, "families-status-summary-yield")
        == "Solved 1 of 1 reviewed (100%) · 2 analysable"
    )
    with user.client:
        families["on_filter"]({"filters": {"Status": ["Solved"]}})
    await wait_for(lambda: [r["FID"] for r in table.rows] == ["1001"])
    await wait_for(
        lambda: (
            _text(user, "families-status-summary-yield")
            == "Solved 1 of 1 reviewed (100%) · 1 analysable"
        )
    )


async def test_the_home_card_shows_the_statuses(user, data_dir):
    datagen.add_family_statuses(data_dir)
    await login(user, "reader")
    await user.should_see(marker="home-status-DEMO-yield", retries=20)
    assert _text(user, "home-status-DEMO-yield") == "Solved 1/1 reviewed (100%)"


async def test_the_report_states_the_status(user, data_dir):
    datagen.add_family_statuses(data_dir)
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=20)
    user.find(marker="family-report-view").click()
    await user.should_see(marker="family-report-frame", retries=60)
    html = user.find(marker="family-report-frame").elements.pop().props["srcdoc"]
    assert "Review status: <strong>Solved (curator, 2026-10-05)</strong>" in html
    # The comment is an internal note.
    assert "SHANK3 de novo" not in html


async def test_the_standalone_family_page_has_the_chip(user, data_dir):
    datagen.add_family_statuses(data_dir)
    await login(user, "curator")
    await user.open("/family/1002")
    await user.should_see(marker="family-status", retries=20)
    assert _text(user, "family-status") == "Status: Not analysable"
    # The curator's tooltip has the comment.
    await user.should_see("tube mix-up at the centre")
