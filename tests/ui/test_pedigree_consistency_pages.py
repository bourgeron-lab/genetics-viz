"""Pedigrees across cohorts: family banner, Families badge, home warning, list, editor.

``datagen.add_cohort_disagreements``: NOHEAD lists 1001 with S3's parents at 0
(DEMO links the trio); cohort ALIAS knows AU-X-1 as AU-77 and S3 alone as K3.
"""

import asyncio
from pathlib import Path

import datagen
import pytest
from ui_helpers import capture_tables, login, open_tab

from genetics_viz.utils import pedigree_consistency as pc

SLOW = 200


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    """The synthetic data with the disagreements, before the app boots."""
    data = datagen.build_data_dir(tmp_path)
    datagen.add_cohort_disagreements(data)
    return data


async def wait_for(condition, tries: int = 200) -> None:
    for _ in range(tries):
        if condition():
            return
        await asyncio.sleep(0.05)
    assert condition()


async def test_a_reader_sees_where_else_the_family_is(user, data_dir):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see(marker="family-consistency", retries=SLOW)
    await user.should_see("Rows differ between cohorts")
    await user.should_see("1001 in NOHEAD")
    await user.should_see("K3 in ALIAS")
    # Readers get no details and no button.
    await user.should_not_see(
        "S3: father S1 (DEMO) / 0 (NOHEAD); mother S2 (DEMO) / 0 (NOHEAD)"
    )
    await user.should_not_see(marker="consistency-explain")


async def test_a_curator_explains_an_alias(user, data_dir):
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/AU-X-1")
    await user.should_see(marker="consistency-explain", retries=SLOW)
    await user.should_see("AU-77 in ALIAS")
    user.find(marker="consistency-explain").click()
    await user.should_see(marker="consistency-dialog", retries=20)
    user.find(marker="consistency-note").type("the centre's own code")
    user.find(marker="consistency-save").click()
    await user.should_see("Explained: the centre's own code", retries=SLOW)
    (ack,) = pc.load_acks(data_dir).values()
    assert (ack.user, ack.note) == ("curator", "the centre's own code")


async def test_the_families_tab_badges_and_the_status_reason(
    user, data_dir, monkeypatch
):
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=20)
    families = next(t for t in tables if "Structure" in t["columns"])
    table = families["table"]
    await wait_for(lambda: any(r.get("Pedigree_badges") for r in table.rows))
    rows = {r["FID"]: r for r in table.rows}
    labels = {
        fid: [b["label"] for b in r["Pedigree_badges"]] for fid, r in rows.items()
    }
    assert "Differs: NOHEAD" in labels["1001"] and "Alias: K3" in labels["1001"]
    assert labels["AU-X-1"] == ["Alias: AU-77"]
    # A high disagreement asks for a check of the status; an alias does not.
    assert "Needs re-review" in rows["1001"]["_status_cats"]
    assert "Needs re-review" not in rows["AU-X-1"]["_status_cats"]


async def test_home_warning_and_the_list(user, data_dir):
    await login(user, "curator")
    await user.should_see(marker="home-consistency", retries=SLOW)
    await user.should_see("3 findings (1 high) over 4 families")
    await user.open("/pedigree-consistency")
    await user.should_see(marker="consistency-count", retries=SLOW)
    await user.should_see("3 to review, 0 explained or informative")


async def test_readers_get_no_home_warning_nor_list(user, data_dir):
    await login(user, "reader")
    await asyncio.sleep(0.3)
    await user.should_not_see(marker="home-consistency")
    await user.open("/pedigree-consistency")
    await user.should_see("Curators and administrators only.")


async def test_the_editor_makes_the_cohorts_agree(user, data_dir):
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=SLOW)
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see(marker="pedigree-edit", retries=SLOW)
    user.find(marker="pedigree-edit").click()
    await user.should_see(marker="ped-reference", retries=SLOW)
    user.find(marker="ped-reference-DEMO").click()
    await user.should_see(marker="ped-cohorts-differ", retries=SLOW)
    await user.should_see("NOHEAD: S3: father 0 -> S1, mother 0 -> S2")
    user.find(marker="ped-comment").type("link the trio in NOHEAD")
    await wait_for(lambda: user.find(marker="ped-save").elements.pop().enabled)
    user.find(marker="ped-save").click()
    nohead = data_dir / "cohorts" / "NOHEAD" / "NOHEAD.pedigree.tsv"
    await wait_for(lambda: "1001\tS3\tS1\tS2" in nohead.read_text())
