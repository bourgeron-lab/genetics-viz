"""The pedigree editor: a curator's dialog from the family's Sample QC tab.

1001 is DEMO's trio (S1, S2, S3), listed in DEMO and NOHEAD. The edit made
here is the one Sample QC suggests for its ``NOT_PARENT`` flag: S3's father
becomes a dummy.
"""

import asyncio
from datetime import datetime
from pathlib import Path

from ui_helpers import capture_tables, login, open_tab

from genetics_viz.utils import pedigree_edit as pe

DEMO = Path("cohorts") / "DEMO" / "DEMO.pedigree.tsv"


async def wait_for(condition, tries: int = 100) -> None:
    for _ in range(tries):
        if condition():
            return
        await asyncio.sleep(0.05)
    assert condition()


def _one(user, marker: str):
    (element,) = user.find(marker=marker).elements
    return element


#: Lazy tabs load off-thread: generous waits, for a busy CI runner.
SLOW = 200


async def open_qc(user, username: str, fid: str = "1001") -> None:
    await login(user, username)
    await user.open(f"/cohort/DEMO/family/{fid}")
    await user.should_see(f"Family: {fid}", retries=SLOW)
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=SLOW)


async def open_editor(user) -> None:
    user.find(marker="pedigree-edit").click()
    await user.should_see(marker="pedigree-editor", retries=SLOW)


async def test_readers_get_no_editor(user):
    await open_qc(user, "reader")
    await user.should_not_see(marker="pedigree-edit")


async def test_a_curator_gives_a_child_a_dummy_father_and_saves(user, data_dir):
    await open_qc(user, "curator")
    await open_editor(user)
    now = _one(user, "ped-edit-now").content
    assert all(f'data-id="{i}"' in now for i in ("S1", "S2", "S3"))

    _one(user, "ped-father-S3").set_value("__new_dummy__")
    await user.should_see("S3: father S1 -> DUMMY_S3_F", retries=10)
    await user.should_see("DUMMY_S3_F: added")
    after = _one(user, "ped-edit-new").content
    assert 'data-id="DUMMY_S3_F"' in after and "node changed" in after
    assert "node changed" not in _one(user, "ped-edit-now").content
    await user.should_see("Rerun ghfc-ngs as usual")
    # A comment is required.
    assert not _one(user, "ped-save").enabled
    user.find(marker="ped-comment").type("S1 is not S3's father (Sample QC)")
    await wait_for(lambda: _one(user, "ped-save").enabled)
    user.find(marker="ped-save").click()

    await wait_for(lambda: "DUMMY_S3_F" in (data_dir / DEMO).read_text())
    assert "1001\tS3\tDUMMY_S3_F\tS2\t1\t2" in (data_dir / DEMO).read_text()
    assert "DUMMY_S3_F" in (data_dir / "cohorts/NOHEAD/NOHEAD.pedigree.tsv").read_text()
    (record,) = pe.read_family_edits(data_dir, "1001")
    assert (record.user, record.comment) == (
        "curator",
        "S1 is not S3's father (Sample QC)",
    )
    assert record.cohorts == ("DEMO", "NOHEAD") and record.rerun == pe.RERUN_NORMAL

    # The page reloads on the new pedigree: the member card has the dummy, and
    # the Sample QC tab says which results still read the older one.
    await user.should_see("Family: 1001", retries=SLOW)
    await user.should_see("DUMMY_S3_F", retries=SLOW)
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=SLOW)
    await user.should_see(marker="pedigree-pending", retries=SLOW)
    await user.should_see("Computed on an older pedigree: Sample QC")
    await user.should_see("Pedigree last edited on")


async def test_errors_block_the_save_and_undo_takes_them_back(user):
    await open_qc(user, "curator")
    await open_editor(user)
    user.find(marker="ped-comment").type("why")
    _one(user, "ped-sex-S2").set_value("1")
    await user.should_see("S2 is S3's mother but is recorded as male", retries=10)
    assert not _one(user, "ped-save").enabled
    user.find(marker="ped-undo").click()
    await user.should_not_see("recorded as male", retries=10)
    await user.should_see("None yet.")
    assert not _one(user, "ped-save").enabled  # nothing changed


async def test_full_sibs_share_dummy_parents(user, data_dir):
    # AU-X-1: S5 (no parents) and S6 (father S5, no mother). Make S5 and S6
    # lack the same parents first: S6's father unknown.
    await open_qc(user, "curator", "AU-X-1")
    await open_editor(user)
    _one(user, "ped-father-S6").set_value("0")
    await user.should_see("S6: father S5 -> 0", retries=10)
    _one(user, "ped-shared-children").set_value(["S5", "S6"])
    _one(user, "ped-shared-role").set_value("father")
    user.find(marker="ped-shared-add").click()
    await user.should_see("S5: father 0 -> DUMMY_S5_F", retries=10)
    await user.should_see("S6: father S5 -> DUMMY_S5_F")
    await user.should_see("DUMMY_S5_F: added")


async def test_the_history_restores_the_rows_before_an_edit(
    user, data_dir, monkeypatch
):
    monkeypatch.setattr(pe, "can_write", lambda: True)
    files = [data_dir / DEMO, data_dir / "cohorts/NOHEAD/NOHEAD.pedigree.tsv"]
    source = pe.load_family(data_dir, "1001", files)
    draft = source.draft()
    draft.set_parent("S3", "father", None)
    draft.add_dummy_parents(["S3"], roles=["father"])
    pe.save_pedigree_edit(source, draft, user="alice", comment="non-paternity")

    await open_qc(user, "curator")
    await user.should_see("Pedigree last edited on", retries=SLOW)
    await open_editor(user)
    _one(user, "ped-history").set_value(True)
    await user.should_see("non-paternity", retries=10)
    (edit_id,) = [r.edit_id for r in pe.read_family_edits(data_dir, "1001")]
    user.find(marker=f"ped-show-{edit_id}").click()
    await user.should_see(marker=f"ped-restore-{edit_id}", retries=30)
    user.find(marker=f"ped-restore-{edit_id}").click()
    await user.should_see("S3: father DUMMY_S3_F -> S1", retries=10)
    await user.should_see("DUMMY_S3_F: removed")


async def test_the_families_tab_shows_edited_families(user, monkeypatch, data_dir):
    monkeypatch.setattr(pe, "can_write", lambda: True)
    source = pe.load_family(data_dir, "1002", [data_dir / DEMO])
    draft = source.draft()
    draft.set_phenotype("S4", "2")
    pe.save_pedigree_edit(
        source, draft, user="alice", comment="affected", now=datetime(2026, 9, 1)
    )
    source = pe.load_family(data_dir, "AU-X-1", [data_dir / DEMO])
    draft = source.draft()
    draft.set_sex("S6", "1")
    pe.save_pedigree_edit(
        source, draft, user="alice", comment="XY", now=datetime(2026, 11, 18)
    )
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    families = next(t for t in tables if "Structure" in t["columns"])
    assert "Pedigree" in families["columns"]
    table = families["table"]
    await wait_for(lambda: any(r.get("Pedigree") for r in table.rows))
    rows = {r["FID"]: r for r in table.rows}
    assert rows["1001"]["Pedigree"] == ""
    # 1002's phenotype edit reruns nothing; AU-X-1's sex edit comes after the
    # last run (2026-10-01): awaiting it.
    assert rows["1002"]["_ped_cats"] == ["Edited"]
    assert rows["AU-X-1"]["_ped_cats"] == ["Edited", "Awaiting a rerun"]
    (badge,) = rows["AU-X-1"]["Pedigree_badges"]
    assert (
        badge["label"] == "2026-11-18" and "No ghfc-ngs run since" in badge["tooltip"]
    )
    with user.client:
        families["on_filter"]({"filters": {"Pedigree": True}})
    await wait_for(lambda: [r["FID"] for r in table.rows] == ["AU-X-1"])


async def test_readers_get_no_pedigree_column(user, monkeypatch):
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    families = next(t for t in tables if "Structure" in t["columns"])
    assert "Pedigree" not in families["columns"]


async def test_sample_qc_suggests_a_dummy_father_and_the_edit_records_the_flag(
    user, data_dir
):
    """1001: S1 is observed unrelated to S3, its declared son (NOT_PARENT)."""
    await open_qc(user, "curator")
    await open_editor(user)
    await user.should_see("Suggestions from Sample QC", retries=10)
    await user.should_see("Sample QC against this pedigree: 1 contradiction now.")
    await user.should_see("S1 – S3: expected parent-offspring, observed unrelated")
    await user.should_see("contradictions 1 → 0")
    user.find(marker="ped-apply-notparent:S3:father").click()
    await user.should_see("S3: father S1 -> DUMMY_S3_F", retries=10)
    await user.should_see(
        "Sample QC against this pedigree: 1 contradiction now; 0 contradictions,"
        " 1 pair related beyond it with your changes."
    )
    user.find(marker="ped-comment").type("Non-paternity, confirmed by the centre")
    await wait_for(lambda: _one(user, "ped-save").enabled)
    user.find(marker="ped-save").click()
    await wait_for(lambda: pe.read_family_edits(data_dir, "1001") != [])
    (record,) = pe.read_family_edits(data_dir, "1001")
    assert record.flags == ("NOT_PARENT",)


async def test_sample_qc_suggests_the_karyotypes_sex(user):
    """AU-X-1: S6 is recorded female, its karyotype is XY (SEX_MISMATCH)."""
    await open_qc(user, "curator", "AU-X-1")
    await open_editor(user)
    await user.should_see("Record S6 as male", retries=10)
    await user.should_see("S6: karyotype XY")
    user.find(marker="ped-apply-sex:S6").click()
    await user.should_see("S6: sex 2 -> 1", retries=10)


async def test_a_swap_with_another_family_is_saved_for_both(user, data_dir):
    """AU-X-1 and 1002 are both in DEMO only: one edit can hold both."""
    await open_qc(user, "curator", "AU-X-1")
    await open_editor(user)
    _one(user, "ped-other-family").set_value("1002")
    user.find(marker="ped-other-open").click()
    await user.should_see("Edit the pedigrees of families AU-X-1 and 1002", retries=30)
    user.find(marker="ped-swap-S6-S4").click()
    await user.should_see("S4: added", retries=10)
    await user.should_see("S6: added")
    await user.should_see(marker="ped-rerun-1002")
    await user.should_see("ghfc-ngs needs --clean-stale-families")
    user.find(marker="ped-comment").type("S6 and S4 were swapped at the lab")
    await wait_for(lambda: _one(user, "ped-save").enabled)
    user.find(marker="ped-save").click()
    await wait_for(lambda: len(pe.read_edit_index(data_dir)) == 2)
    text = (data_dir / DEMO).read_text()
    assert "AU-X-1\tS4\tS5\t0\t2\t1" in text and "1002\tS6\t0\t0\t2.0\t-9.0" in text
    assert {r.rerun for r in pe.read_edit_index(data_dir)} == {pe.RERUN_CLEAN}

    # The cohort's Workflow tab lists both, for the clean. A fresh cached
    # parameters reference: the tab reads it instead of fetching GitHub.
    import sys

    docs = sys.modules["genetics_viz.utils.ghfc_ngs_docs"]
    docs._cached = docs.FetchedDoc(
        markdown="# Stub reference", error=None, fetched_at=datetime.now()
    )
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-tabs", "Workflow")
    # The tab warms its caches off-thread first: slow on a busy CI runner.
    await user.should_see("Pipeline status", retries=SLOW)
    await user.should_see(marker="workflow-pedigree-waiting", retries=100)
    await user.should_see("--clean-stale-families for 1002, AU-X-1")


def add_cross_pair(data_dir: Path, fa, a, fb, b, observed: str, phi: float) -> None:
    """One more cross-family pair in DEMO's cohort Sample QC, as the pipeline writes it."""
    import datagen

    path = data_dir / "cohorts" / "DEMO" / "qc" / f"DEMO.{datagen.QC_LABEL}.pairs.tsv"
    row = datagen._qc_pair_row(
        fa, a, fb, b, phi, 0.2, 0.0, "unrelated", observed, None, "mismatch"
    )
    with open(path, "a") as f:
        f.write("\t".join(map(str, row)) + "\n")


async def test_merge_candidates_on_the_cohort_sample_qc_tab(user, data_dir):
    # 1001 and AU-X-1 (S2–S5, second degree) are in different cohort pedigrees;
    # 1002 and AU-X-1, both in DEMO only, are made full sibs here.
    add_cross_pair(data_dir, "1002", "S4", "AU-X-1", "S6", "full-sibs", 0.25)
    await login(user, "curator")
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-tabs", "Sample QC")
    await user.should_see(marker="merge-candidates", retries=30)
    await user.should_see("1001 + 1002 + AU-X-1")
    await user.should_see("S4 (1002) – S6 (AU-X-1): observed full-sibs")
    button = _one(user, "merge-candidate-1001+1002+AU-X-1")
    assert not button.enabled  # 1001 is also in NOHEAD


async def test_readers_get_no_merge_candidates(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-tabs", "Sample QC")
    await user.should_see("Relatives across families", retries=30)
    await user.should_not_see(marker="merge-candidates")


async def test_two_families_merged_under_a_new_id(user, data_dir):
    await open_qc(user, "curator", "AU-X-1")
    await open_editor(user)
    _one(user, "ped-merge-with").set_value(["1002"])
    user.find(marker="ped-merge-start").click()
    await user.should_see("Merge families AU-X-1, 1002", retries=30)
    assert _one(user, "ped-merge-fid").value == "AU-X-1_1002"
    await user.should_see("AU-X-1, 1002 merged into AU-X-1_1002 (a new family)")
    await user.should_see(marker="ped-merge-curation")  # S6's acknowledgement
    _one(user, "ped-merge-fid").set_value("1001")
    await user.should_see(
        "Family 1001 already exists in the cohort pedigree", retries=10
    )
    _one(user, "ped-merge-fid").set_value("FAMX")
    await user.should_see("AU-X-1, 1002 merged into FAMX (a new family)", retries=10)
    user.find(marker="ped-comment").type("One family entered under two IDs")
    await wait_for(lambda: _one(user, "ped-save").enabled)
    user.find(marker="ped-save").click()

    # The merged family's page opens, on the new pedigree.
    await user.should_see("Family: FAMX", retries=SLOW)
    text = (data_dir / DEMO).read_text()
    assert "FAMX\tS4" in text and "FAMX\tS6" in text and "AU-X-1\t" not in text
    acks = (data_dir / "validations" / "sampleqc.tsv").read_text()
    assert "SEX_MISMATCH\tFAMX\tS6" in acks

    # The Workflow tab names the old directories, and FAMX as waiting.
    import sys

    docs = sys.modules["genetics_viz.utils.ghfc_ngs_docs"]
    docs._cached = docs.FetchedDoc(
        markdown="# Stub reference", error=None, fetched_at=datetime.now()
    )
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-tabs", "Workflow")
    await user.should_see("Pipeline status", retries=SLOW)
    await user.should_see(marker="workflow-orphaned-families", retries=100)
    await user.should_see("AU-X-1 → FAMX", retries=10)
    await user.should_see(marker="workflow-pedigree-waiting")
