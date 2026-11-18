"""Pedigree edits: the draft, its checks, the files written, the history.

The data are ``datagen``'s DEMO cohort: 1001 a trio (S1, S2, S3), 1002 one
woman written ``2.0``/``-9.0``, AU-X-1 a father and daughter whose missing
parents are written ``-9``, empty and ``0``. The cohort pedigree has a header.
"""

import os
from datetime import datetime
from pathlib import Path

import datagen
import pytest

from genetics_viz.models import DataStore
from genetics_viz.utils import pedigree_edit as pe
from genetics_viz.utils.sharding import get_family_path

NOW = datetime(2026, 11, 18, 9, 30, 0)


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    return datagen.build_data_dir(tmp_path)


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(pe, "can_write", lambda: True)


def cohort_file(data_dir: Path, name: str = "DEMO") -> Path:
    return data_dir / "cohorts" / name / f"{name}.pedigree.tsv"


def source(data_dir: Path, fid: str = "1001", *names: str) -> pe.FamilySource:
    files = [cohort_file(data_dir, n) for n in (names or ("DEMO",))]
    return pe.load_family(data_dir, fid, files)


def draft_of(*rows) -> pe.PedigreeDraft:
    return pe.PedigreeDraft("F", [pe.Member(*r) for r in rows])


def messages(found, level=None):
    return [p.message for p in found if level is None or p.level == level]


# --- Codes, as ghfc-ngs reads them -------------------------------------------------


@pytest.mark.parametrize(
    "cells, row",
    [
        (["F", "K", "-9", "", "1.0", "2"], ("F", "K", "0", "0", "1")),
        (["F", "K", ".", "0", "2", "-9"], ("F", "K", "0", "0", "2")),
        (["F", "K", "D", "M", "x", "1"], ("F", "K", "D", "M", "0")),
        # The app reads 0.0 as missing; the pipeline reads it as an ID.
        (["F", "K", "0.0", "0", "2.00", "1"], ("F", "K", "0.0", "0", "2")),
    ],
)
def test_snapshot_row_is_the_pipelines(cells, row):
    assert pe.snapshot_row(cells) == row


# --- Operations --------------------------------------------------------------------


def test_full_sibs_share_one_dummy_parent_half_sibs_one_each():
    draft = draft_of(
        ("M", "0", "0", "2", "1"),
        ("K1", "0", "M", "1", "2"),
        ("K2", "0", "M", "2", "1"),
        ("K3", "0", "M", "2", "1"),
    )
    assert draft.add_dummy_parents(["K1", "K2"], roles=["father"]) == ["DUMMY_K1_F"]
    assert draft.add_dummy_parents(["K3"], roles=["father"]) == ["DUMMY_K3_F"]
    fathers = {
        m.iid: m.parent("father") for m in draft.members if m.iid.startswith("K")
    }
    assert fathers == {"K1": "DUMMY_K1_F", "K2": "DUMMY_K1_F", "K3": "DUMMY_K3_F"}
    dummy = draft.get("DUMMY_K1_F")
    assert (dummy.sex, dummy.father, dummy.mother, dummy.phenotype) == (
        "1",
        "0",
        "0",
        "0",
    )
    # A parent already known is left alone: K1 keeps its mother.
    assert draft.add_dummy_parents(["K1"]) == []


def test_a_dummys_dummy_parent_extends_its_path():
    assert pe.dummy_parent_id("S9", "mother") == "DUMMY_S9_M"
    assert pe.dummy_parent_id("DUMMY_S9_M", "father") == "DUMMY_S9_MF"


def test_swap_exchanges_places_and_parent_references():
    draft = draft_of(
        ("A", "0", "0", "1", "1"),
        ("B", "0", "0", "2", "1"),
        ("C", "A", "B", "1", "2"),
    )
    draft.swap("A", "C")
    rows = {m.iid: (m.father, m.mother, m.sex, m.phenotype) for m in draft.members}
    assert rows == {
        "C": ("0", "0", "1", "1"),
        "B": ("0", "0", "2", "1"),
        "A": ("C", "B", "1", "2"),
    }
    with pytest.raises(pe.PedigreeEditError):
        draft.swap("A", "DUMMY_X_F")


def test_an_implied_dummy_gets_a_row_once_edited():
    draft = draft_of(("K", "DUMMY_K_F", "0", "1", "2"))
    assert [(m.iid, m.sex) for m in draft.implied] == [("DUMMY_K_F", "1")]
    draft.set_parent("DUMMY_K_F", "father", None)  # no change, but now a row
    assert "DUMMY_K_F" in draft and draft.implied == []


def test_a_dummy_parent_cannot_be_removed():
    draft = draft_of(
        ("K", "DUMMY_K_F", "0", "1", "2"), ("DUMMY_K_F", "0", "0", "1", "0")
    )
    with pytest.raises(pe.PedigreeEditError, match="parent of K"):
        draft.remove_dummy("DUMMY_K_F")
    draft.set_parent("K", "father", "0")
    draft.remove_dummy("DUMMY_K_F")
    assert draft.ids == ["K"]


# --- Checks ------------------------------------------------------------------------


def test_errors_block_and_warnings_do_not():
    original = draft_of(
        ("D", "0", "0", "1", "1"),
        ("M", "0", "0", "2", "1"),
        ("K", "D", "M", "1", "2"),
    )
    draft = original.copy()
    draft.set_sex("D", "2")
    found = pe.problems(original, draft)
    assert messages(found, pe.ERROR) == ["D is K's father but is recorded as female"]
    assert pe.has_errors(found)

    draft = original.copy()
    draft.set_parent("D", "father", "K")  # K is D's son and father
    assert any("own ancestor" in m for m in messages(pe.problems(original, draft)))

    draft = original.copy()
    draft.set_parent("K", "mother", "D")
    assert "K: D is both father and mother" in messages(pe.problems(original, draft))


def test_new_dummies_are_named_by_convention_and_unique_in_the_file():
    original = draft_of(("K", "0", "0", "1", "2"))
    draft = original.copy()
    draft.add_dummy("DUMMY_whatever", "1")
    draft.add_dummy("DUMMY_S1_F", "1")
    draft.set_parent("K", "father", "DUMMY_whatever")
    errors = messages(pe.problems(original, draft, other_ids={"DUMMY_S1_F"}), pe.ERROR)
    assert errors == [
        "DUMMY_whatever: a new dummy is named DUMMY_<child>_<F|M…>",
        "DUMMY_S1_F is already a member of another family",
    ]


def test_sequenced_members_cannot_change_yet():
    original = draft_of(("A", "0", "0", "1", "1"), ("B", "0", "0", "2", "1"))
    draft = pe.PedigreeDraft("F", [original.members[0]])
    (error,) = messages(pe.problems(original, draft), pe.ERROR)
    assert "removed B" in error and "another family open alongside" in error
    assert pe.rerun_needed(original, draft) == pe.RERUN_CLEAN


def test_warnings_name_undecided_sibships_and_unreadable_codes():
    original = draft_of(
        ("D", "0", "0", "1", "1"),
        ("K1", "D", "0", "1", "2"),
        ("K2", "D", "-9", "2", "1"),
        ("K3", "0.0", "0", "0", "1"),
        ("DUMMY_K1_M", "0", "0", "2", "0"),
    )
    warnings = messages(pe.problems(original, original), pe.WARNING)
    assert (
        "K3's father is written 0.0, which ghfc-ngs reads as an ID, not as unknown: set it to unknown"
        in warnings
    )
    assert (
        "K3's father 0.0 has no row: the pipeline reads them as a founder" in warnings
    )
    assert "DUMMY_K1_M is no one's parent" in warnings
    assert "K3's sex is unknown" in warnings
    assert any(
        w.startswith("K1, K2 share their father D and have no mother") for w in warnings
    )
    assert not pe.has_errors(pe.problems(original, original))


# --- Changes and rerun -------------------------------------------------------------


def test_changes_are_worded_as_the_pipelines_run_summary():
    original = draft_of(("S1", "0", "0", "1.0", "1"), ("S3", "S1", "0", "1", "2"))
    draft = original.copy()
    draft.add_dummy_parents(["S3"], roles=["mother"])
    draft.set_phenotype("S1", "2")
    draft.set_sex("S1", "1")  # 1.0 -> 1 is no change
    assert pe.describe_changes(original, draft) == [
        "DUMMY_S3_M: added",
        "S1: phenotype 1 -> 2",
        "S3: mother 0 -> DUMMY_S3_M",
    ]
    assert pe.rerun_needed(original, draft) == pe.RERUN_NORMAL


def test_a_phenotype_edit_reruns_nothing():
    original = draft_of(("S1", "0", "0", "1", "1"))
    draft = original.copy()
    draft.set_phenotype("S1", "2")
    assert pe.rerun_needed(original, draft) == pe.RERUN_NONE


# --- The cohort file -----------------------------------------------------------------


def test_only_the_familys_lines_change(data_dir):
    path = cohort_file(data_dir)
    path.write_text(
        "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\n"
        "# a comment\n"
        "1002\tS4\t0\t0\t2.0\t-9.0\n"
        "1001\tS1\t0\t0\t1\t1\n"
        "1001\tS2\t0\t0\t2\t1\n"
        "\n"
        "1001\tS3\tS1\t0\t1\t2\textra\n"
        "1001\tS3\tS1\t0\t1\t2\textra\n"
        "AU-X-1\tS5\t-9\t\t1\t2\n"
        "AU-X-1\tS6\tS5\t0\t2\t1"  # no final newline
    )
    src = source(data_dir)
    assert src.original.ids == ["S1", "S2", "S3"]
    draft = src.draft()
    draft.set_parent("S3", "mother", "S2")
    draft.add_dummy_parents(["S1"], roles=["father"])
    text = pe.render_cohort_text(src.files[0], src.original, draft)
    assert text == (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\n"
        "# a comment\n"
        "1002\tS4\t0\t0\t2.0\t-9.0\n"
        "1001\tS1\tDUMMY_S1_F\t0\t1\t1\n"
        "1001\tS2\t0\t0\t2\t1\n"
        "\n"
        "1001\tS3\tS1\tS2\t1\t2\textra\n"  # extra kept, the repeat dropped
        "1001\tDUMMY_S1_F\t0\t0\t1\t0\n"  # after the family's last line
        "AU-X-1\tS5\t-9\t\t1\t2\n"
        "AU-X-1\tS6\tS5\t0\t2\t1"
    )
    assert pe.family_file_text(text, "1001") == (
        "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\n"
        "1001\tS1\tDUMMY_S1_F\t0\t1\t1\n"
        "1001\tS2\t0\t0\t2\t1\n"
        "1001\tS3\tS1\tS2\t1\t2\textra\n"
        "1001\tDUMMY_S1_F\t0\t0\t1\t0\n"
    )


def test_a_rewritten_row_always_has_a_phenotype(data_dir):
    path = cohort_file(data_dir, "NOHEAD")
    path.write_text("AU-X-1\tS5\t-9\t\t1\t\nAU-X-1\tS6\tS5\t0\t2\t1\n")
    src = source(data_dir, "AU-X-1", "NOHEAD")
    draft = src.draft()
    draft.set_sex("S5", "1")
    draft.set_parent("S5", "mother", "0")
    text = pe.render_cohort_text(src.files[0], src.original, draft)
    assert text.splitlines()[0] == "AU-X-1\tS5\t-9\t0\t1\t0"
    # No header in the cohort file, none in the family file.
    assert pe.family_file_text(text, "AU-X-1").startswith("AU-X-1\tS5")


def test_cohorts_that_disagree_are_made_to_agree(data_dir, curator):
    """NOHEAD leaves S3's parents at 0 (PMS's C0733-011-021): editing from
    DEMO's rows writes them into NOHEAD too, NOHEAD's phenotypes kept."""
    nohead = cohort_file(data_dir, "NOHEAD")  # DEMO's first three rows: 1001
    assert not source(data_dir, "1001", "DEMO", "NOHEAD").differs
    nohead.write_text(
        nohead.read_text()
        .replace("S1\tS2\t1\t2", "0\t0\t1\t2")
        .replace("S2\t0\t0\t2\t1", "S2\t0\t0\t2\t2")
    )
    src = pe.load_family(
        data_dir,
        "1001",
        [cohort_file(data_dir, "NOHEAD"), cohort_file(data_dir, "DEMO")],
        reference=cohort_file(data_dir, "DEMO"),
    )
    assert src.differs and src.reference_cohort == "DEMO"
    assert src.file_changes(src.draft()) == [
        "NOHEAD: S3: father 0 -> S1, mother 0 -> S2"
    ]
    assert src.rerun(src.draft()) == pe.RERUN_NORMAL
    (record,) = pe.save_cross_family_edit(
        [(src, src.draft())], user="alice", comment="trio links", now=NOW
    )
    assert record.changes == ("NOHEAD: S3: father 0 -> S1, mother 0 -> S2",)
    text = nohead.read_text()
    assert "1001\tS3\tS1\tS2\t1\t2" in text
    # NOHEAD says S2 is affected: a cohort's own phenotype stays.
    assert "1001\tS2\t0\t0\t2\t2" in text
    assert "1001\tS2\t0\t0\t2\t1" in cohort_file(data_dir, "DEMO").read_text()


def test_a_phenotype_edit_reaches_every_cohort(data_dir, curator):
    src = source(data_dir, "1001", "DEMO", "NOHEAD")
    draft = src.draft()
    draft.set_phenotype("S1", "2")
    pe.save_pedigree_edit(src, draft, user="alice", comment="affected", now=NOW)
    for name in ("DEMO", "NOHEAD"):
        assert "1001\tS1\t0\t0\t1\t2" in cohort_file(data_dir, name).read_text()


def test_the_cohorts_holding_a_family(data_dir):
    store = DataStore(data_dir)
    store.load()
    files = pe.family_cohort_files(store, "1001")
    assert sorted(p.parent.name for p in files) == ["DEMO", "NOHEAD"]


# --- Saving --------------------------------------------------------------------------


def edit_1001(data_dir: Path) -> tuple:
    src = source(data_dir, "1001", "DEMO", "NOHEAD")
    draft = src.draft()
    draft.set_parent("S3", "father", None)
    draft.add_dummy_parents(["S3"], roles=["father"])
    return src, draft


def others(text: str) -> list:
    return [line for line in text.splitlines() if not line.startswith("1001")]


def test_a_save_writes_both_cohorts_the_family_file_and_the_history(data_dir, curator):
    demo_before = cohort_file(data_dir).read_text()
    src, draft = edit_1001(data_dir)
    record = pe.save_pedigree_edit(
        src,
        draft,
        user="alice",
        comment="S1 is not S3's father (QC)",
        flags=["NOT_PARENT"],
        now=NOW,
    )
    assert record.edit_id == "20261118T093000"
    assert record.changes == ("DUMMY_S3_F: added", "S3: father S1 -> DUMMY_S3_F")
    assert record.cohorts == ("DEMO", "NOHEAD") and record.rerun == pe.RERUN_NORMAL

    demo = cohort_file(data_dir).read_text()
    assert "1001\tS3\tDUMMY_S3_F\tS2\t1\t2\n1001\tDUMMY_S3_F\t0\t0\t1\t0\n" in demo
    # Other families untouched, line for line.
    assert others(demo) == others(demo_before)
    assert "DUMMY_S3_F" in cohort_file(data_dir, "NOHEAD").read_text()

    family_file = get_family_path(data_dir, "1001") / "1001.pedigree.tsv"
    assert family_file.read_text() == pe.family_file_text(demo, "1001")
    assert (
        pe.edit_rows(data_dir, "1001", record.edit_id, "after")
        == family_file.read_text()
    )
    assert "1001\tS3\tS1\tS2" in pe.edit_rows(
        data_dir, "1001", record.edit_id, "before"
    )

    (logged,) = pe.read_family_edits(data_dir, "1001")
    assert logged == record
    assert pe.last_edits(data_dir)["1001"].comment == "S1 is not S3's father (QC)"
    # No temporary file left beside a target.
    assert not list(cohort_file(data_dir).parent.glob(".DEMO.pedigree.tsv.*"))

    # The app reads the new pedigree.
    store = DataStore(data_dir)
    store.load()
    assert (
        store.cohorts["DEMO"].families["1001"].get_sample("S3").father_id
        == "DUMMY_S3_F"
    )


def test_a_second_edit_the_same_second_gets_its_own_id(data_dir, curator):
    src, draft = edit_1001(data_dir)
    pe.save_pedigree_edit(src, draft, user="alice", comment="first", now=NOW)
    src = source(data_dir, "1001", "DEMO", "NOHEAD")
    draft = src.draft()
    draft.set_phenotype("S1", "2")
    record = pe.save_pedigree_edit(src, draft, user="bob", comment="second", now=NOW)
    assert record.edit_id == "20261118T093000-2" and record.rerun == pe.RERUN_NONE
    assert [r.user for r in pe.read_family_edits(data_dir, "1001")] == ["alice", "bob"]
    assert pe.last_edits(data_dir)["1001"].user == "bob"


def nothing_written(data_dir: Path, before: str) -> bool:
    return (
        cohort_file(data_dir).read_text() == before
        and not (get_family_path(data_dir, "1001") / pe.EDITS_DIRNAME).exists()
        and not pe.index_file(data_dir).exists()
    )


def test_a_viewer_cannot_save(data_dir, monkeypatch):
    monkeypatch.setattr(pe, "can_write", lambda: False)
    before = cohort_file(data_dir).read_text()
    src, draft = edit_1001(data_dir)
    with pytest.raises(PermissionError):
        pe.save_pedigree_edit(src, draft, user="eve", comment="x", now=NOW)
    assert nothing_written(data_dir, before)


def write_progress(cohort_dir: Path, status: str, minutes_ago: float) -> None:
    """``.ghfc-ngs.progress.json`` as ghfc-ngs's RunProgress writes it."""
    import json
    from datetime import timedelta

    updated = datetime.now().astimezone() - timedelta(minutes=minutes_ago)
    (cohort_dir / pe.PROGRESS_FILENAME).write_text(
        json.dumps(
            {
                "schema_version": 1,
                "run_name": "tender_allen",
                "slurm_job_id": "23792765",
                "status": status,
                "updated_at": updated.isoformat(timespec="seconds"),
            }
        )
    )


@pytest.mark.parametrize(
    "status, minutes_ago, live",
    [
        ("running", 1, True),
        # The file stays after a run: its last word is success or failed.
        ("success", 1, False),
        ("failed", 1, False),
        # A killed run leaves "running" with an updated_at that stops moving.
        ("running", 60, False),
    ],
)
def test_a_run_is_live_only_while_its_progress_file_moves(
    tmp_path, status, minutes_ago, live
):
    write_progress(tmp_path, status, minutes_ago)
    found = pe.live_run(tmp_path)
    if live:
        assert found == "run tender_allen, Slurm job 23792765, updated 1 min ago"
    else:
        assert found is None


def test_no_or_an_unreadable_progress_file_is_no_run(tmp_path):
    assert pe.live_run(tmp_path) is None
    (tmp_path / pe.PROGRESS_FILENAME).write_text("{not json")
    assert pe.live_run(tmp_path) is None


def test_a_finished_runs_progress_file_does_not_block_a_save(data_dir, curator):
    write_progress(cohort_file(data_dir).parent, "success", minutes_ago=600)
    src, draft = edit_1001(data_dir)
    record = pe.save_pedigree_edit(src, draft, user="alice", comment="why", now=NOW)
    assert record.rerun == pe.RERUN_NORMAL


@pytest.mark.parametrize(
    "case, match",
    [
        ("no comment", "comment is required"),
        ("no change", "Nothing changed"),
        ("error", "recorded as male"),
        ("running", "in progress"),
        ("changed", "changed since the editor was opened"),
    ],
)
def test_a_refused_save_writes_nothing(data_dir, curator, case, match):
    before = cohort_file(data_dir).read_text()
    src, draft = edit_1001(data_dir)
    comment = "why"
    if case == "no comment":
        comment = "  "
    elif case == "no change":
        draft = src.draft()
    elif case == "error":
        draft.set_sex("S2", "1")
    elif case == "running":
        write_progress(cohort_file(data_dir).parent, "running", minutes_ago=1)
    elif case == "changed":
        cohort_file(data_dir, "NOHEAD").write_text(
            cohort_file(data_dir, "NOHEAD").read_text() + "# touched\n"
        )
    with pytest.raises(pe.PedigreeEditError, match=match):
        pe.save_pedigree_edit(src, draft, user="alice", comment=comment, now=NOW)
    assert nothing_written(data_dir, before)


@pytest.mark.skipif(os.geteuid() == 0, reason="root writes anywhere")
def test_an_unwritable_directory_stops_the_save_before_any_change(data_dir, curator):
    before = cohort_file(data_dir).read_text()
    src, draft = edit_1001(data_dir)
    nohead_dir = cohort_file(data_dir, "NOHEAD").parent
    nohead_dir.chmod(0o555)
    try:
        with pytest.raises(pe.PedigreeEditError, match="Cannot write"):
            pe.save_pedigree_edit(src, draft, user="alice", comment="why", now=NOW)
    finally:
        nohead_dir.chmod(0o755)
    assert cohort_file(data_dir).read_text() == before
    assert not list(cohort_file(data_dir).parent.glob(".DEMO.pedigree.tsv.*"))


def test_the_cohort_files_mode_is_kept(data_dir, curator):
    cohort_file(data_dir).chmod(0o640)
    src, draft = edit_1001(data_dir)
    pe.save_pedigree_edit(src, draft, user="alice", comment="why", now=NOW)
    assert cohort_file(data_dir).stat().st_mode & 0o777 == 0o640


# --- Steps waiting for a rerun ---------------------------------------------------------


def write_snapshot(path: Path, step: str, fid: str, rows: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#inputs\tv1\n#step\t{step}\n#label\tL\n#family\t{fid}\n" + rows)


def test_pending_steps_follow_the_pipelines_snapshots(data_dir, curator):
    datagen.add_denovo_calls(data_dir)
    fam = get_family_path(data_dir, "1001")
    trio = [r for r in datagen.DEMO_PEDIGREE if r[0] == "1001"]
    # Wombat sorts its rows by sample; the STR table leaves dummies out.
    write_snapshot(
        fam / "wombat" / "inputs" / "1001.rare.vep1.annotated.default.inputs.tsv",
        "wombat", "1001", datagen.snapshot_rows(sorted(trio, key=lambda r: r[1], reverse=True)),
    )  # fmt: skip
    write_snapshot(
        fam / "str" / "expansionhunter" / "1001.L.inputs.tsv",
        "expansionhunter", "1001", datagen.snapshot_rows(trio),
    )  # fmt: skip
    text = cohort_file(data_dir).read_text()
    assert pe.pending_steps(data_dir, "1001", text) == []

    # A dummy no one names yet: every snapshot but the STR one holds it.
    src = source(data_dir, "1001", "DEMO", "NOHEAD")
    draft = src.draft()
    draft.add_dummy("DUMMY_S1_F", "1")
    pe.save_pedigree_edit(src, draft, user="alice", comment="why", now=NOW)
    pending = pe.pending_steps(data_dir, "1001", cohort_file(data_dir).read_text())
    assert sorted(p.step for p in pending) == ["denovo", "sampleqc", "wombat"]

    # A sex edit reaches the STR snapshot too; the history files are no step.
    src = source(data_dir, "1001", "DEMO", "NOHEAD")
    draft = src.draft()
    draft.set_sex("S3", "2")
    pe.save_pedigree_edit(src, draft, user="alice", comment="why", now=NOW)
    pending = pe.pending_steps(data_dir, "1001", cohort_file(data_dir).read_text())
    assert sorted(p.step for p in pending) == [
        "denovo",
        "expansionhunter",
        "sampleqc",
        "wombat",
    ]


def test_datagen_snapshots_match_the_pedigree(data_dir):
    """Each family's snapshots were written from DEMO's rows: nothing waits."""
    datagen.add_genomeqc_tables(data_dir)
    text = cohort_file(data_dir).read_text()
    for fid in ("1001", "1002", "AU-X-1"):
        assert pe.pending_steps(data_dir, fid, text) == [], fid


# --- Across families ---------------------------------------------------------------


def test_a_moved_parent_leaves_a_shared_dummy_behind():
    a = pe.PedigreeDraft(
        "A",
        [
            pe.Member("D", "DUMMY_D_F", "0", "1", "1"),
            pe.Member("DUMMY_D_F", "0", "0", "1", "0"),
            pe.Member("M", "0", "0", "2", "1"),
            pe.Member("K1", "D", "M", "1", "2"),
            pe.Member("K2", "D", "M", "2", "1"),
        ],
    )
    b = pe.PedigreeDraft("B", [pe.Member("Y", "0", "0", "2", "1")])
    pe.move_member(a, b, "D")
    assert a.ids == ["M", "K1", "K2", "DUMMY_K1_F"]  # D's own dummy father gone
    assert {a.get(k).father for k in ("K1", "K2")} == {"DUMMY_K1_F"}
    assert b.get("D") == pe.Member("D", "0", "0", "1", "1")
    with pytest.raises(pe.PedigreeEditError):
        pe.move_member(a, b, "DUMMY_K1_F")


def test_a_swap_across_families_exchanges_places():
    a = pe.PedigreeDraft(
        "A", [pe.Member("X", "0", "0", "1", "1"), pe.Member("K", "X", "0", "1", "2")]
    )
    b = pe.PedigreeDraft(
        "B", [pe.Member("P", "0", "0", "1", "1"), pe.Member("Y", "P", "0", "1", "2")]
    )
    pe.swap_across(a, b, "X", "Y")
    assert [(m.iid, m.father) for m in a.members] == [("Y", "0"), ("K", "Y")]
    assert [(m.iid, m.father) for m in b.members] == [("P", "0"), ("X", "P")]


def test_cross_problems_check_the_members_as_a_whole():
    a = pe.PedigreeDraft("A", [pe.Member("X", "0", "0", "1", "1")])
    b = pe.PedigreeDraft("B", [pe.Member("Y", "0", "0", "2", "1")])
    moved_a, moved_b = a.copy(), b.copy()
    pe.move_member(moved_a, moved_b, "X")
    errors = [p.message for p in pe.cross_problems([a, b], [moved_a, moved_b])]
    assert errors == ["A would have no sequenced member"]
    lost = b.copy()
    lost.members = []
    errors = [p.message for p in pe.cross_problems([a, b], [a.copy(), lost])]
    assert "Sequenced members cannot be removed or added here: Y" in errors


def test_a_move_between_families_is_saved_once_for_both(data_dir, curator):
    demo = [cohort_file(data_dir)]
    a = pe.load_family(data_dir, "1001", demo)
    b = pe.load_family(data_dir, "AU-X-1", demo)
    da, db = a.draft(), b.draft()
    pe.move_member(da, db, "S3")
    records = pe.save_cross_family_edit(
        [(a, da), (b, db)], user="alice", comment="S3 belongs to AU-X-1", now=NOW
    )
    assert [(r.family_id, r.rerun) for r in records] == [
        ("1001", pe.RERUN_CLEAN),
        ("AU-X-1", pe.RERUN_CLEAN),
    ]
    assert records[0].changes == ("S3: removed",)
    assert records[1].changes == ("S3: added",)
    assert records[0].comment == "S3 belongs to AU-X-1 (with AU-X-1)"
    assert records[0].edit_id == records[1].edit_id
    text = cohort_file(data_dir).read_text()
    assert "1001\tS3" not in text and "AU-X-1\tS3\t0\t0\t1\t2\n" in text
    store = DataStore(data_dir)
    store.load()
    families = store.cohorts["DEMO"].families
    assert "S3" not in [s.sample_id for s in families["1001"].samples]
    assert "S3" in [s.sample_id for s in families["AU-X-1"].samples]
    family = get_family_path(data_dir, "AU-X-1") / "AU-X-1.pedigree.tsv"
    assert "S3" in family.read_text()
    assert [r.family_id for r in pe.read_edit_index(data_dir)] == ["1001", "AU-X-1"]


def test_families_in_different_cohorts_are_edited_apart(data_dir, curator):
    a = pe.load_family(
        data_dir, "1001", [cohort_file(data_dir), cohort_file(data_dir, "NOHEAD")]
    )
    b = pe.load_family(data_dir, "AU-X-1", [cohort_file(data_dir)])
    da, db = a.draft(), b.draft()
    pe.move_member(da, db, "S3")
    with pytest.raises(pe.PedigreeEditError, match="not in the same cohort pedigrees"):
        pe.save_cross_family_edit(
            [(a, da), (b, db)], user="alice", comment="x", now=NOW
        )


def write_bcf(data_dir: Path, fid: str, samples) -> None:
    """A joint call with only a header: the samples ghfc-ngs called together."""
    import pysam

    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000)
    for s in samples:
        header.add_sample(s)
    path = get_family_path(data_dir, fid) / "vcfs" / f"{fid}.norm.bcf"
    path.parent.mkdir(parents=True, exist_ok=True)
    with pysam.VariantFile(str(path), "wb", header=header):
        pass


def test_waiting_families_and_the_clean_they_need(data_dir, curator):
    write_bcf(data_dir, "1001", ["S1", "S2", "S3"])
    demo = [cohort_file(data_dir)]
    a = pe.load_family(data_dir, "1001", demo)
    b = pe.load_family(data_dir, "AU-X-1", demo)
    da, db = a.draft(), b.draft()
    pe.move_member(da, db, "S3")
    pe.save_cross_family_edit([(a, da), (b, db)], user="alice", comment="move", now=NOW)
    c = pe.load_family(data_dir, "1002", demo)
    dc = c.draft()
    dc.set_sex("S4", "1")
    pe.save_pedigree_edit(c, dc, user="bob", comment="sex", now=NOW)

    store = DataStore(data_dir)
    store.load()
    members = {
        fid: [s.sample_id for s in f.samples]
        for fid, f in store.cohorts["DEMO"].families.items()
    }
    run = datetime(2026, 10, 1, 10, 0).astimezone()
    waiting = {w.family_id: w for w in pe.waiting_families(data_dir, members, run)}
    assert set(waiting) == {"1001", "1002", "AU-X-1"}
    # 1001's joint call still has S3: the drift ghfc-ngs will see.
    assert (waiting["1001"].extra, waiting["1001"].needs_clean) == (("S3",), True)
    # AU-X-1 has no joint call here: the edit's own record says so.
    assert waiting["AU-X-1"].missing is None and waiting["AU-X-1"].needs_clean
    assert waiting["1002"].rerun == pe.RERUN_NORMAL and not waiting["1002"].needs_clean
    # A run started after the edits: nothing waits.
    assert (
        pe.waiting_families(data_dir, members, datetime(2027, 1, 1).astimezone()) == []
    )
