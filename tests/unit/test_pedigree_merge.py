"""Merging families: the pedigree lines, the target ID, the curation re-keyed.

1002 (S4) and AU-X-1 (S5, S6) are both in DEMO only, so they can be merged in
one edit; 1001 is also in NOHEAD. ``validations/sampleqc.tsv`` holds an
acknowledgement for AU-X-1 and a cross-family one for ``1001,AU-X-1``.
"""

from datetime import datetime
from pathlib import Path

import datagen
import pytest

from genetics_viz.models import DataStore
from genetics_viz.utils import pedigree_edit as pe
from genetics_viz.utils import pedigree_merge as pm
from genetics_viz.utils.sharding import get_family_path
from genetics_viz.utils.tsv import append_tsv_row

NOW = datetime(2026, 11, 18, 12, 0, 0)


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    return datagen.build_data_dir(tmp_path)


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(pe, "can_write", lambda: True)


def demo(data_dir: Path) -> Path:
    return data_dir / "cohorts" / "DEMO" / "DEMO.pedigree.tsv"


def merge_of(data_dir: Path, *fids: str) -> pm.MergeSource:
    return pm.load_merge(data_dir, list(fids), [[demo(data_dir)] for _ in fids])


def test_the_suggested_id_joins_the_families():
    assert pm.suggested_fid(["1002", "AU-X-1"]) == "1002_AU-X-1"


def test_families_of_other_cohort_pedigrees_cannot_be_merged(data_dir):
    nohead = data_dir / "cohorts" / "NOHEAD" / "NOHEAD.pedigree.tsv"
    with pytest.raises(pe.PedigreeEditError, match="not in the same cohort"):
        pm.load_merge(
            data_dir, ["1001", "1002"], [[demo(data_dir), nohead], [demo(data_dir)]]
        )
    with pytest.raises(pe.PedigreeEditError, match="two families or more"):
        pm.load_merge(data_dir, ["1002"], [[demo(data_dir)]])


@pytest.mark.parametrize(
    "target, error",
    [
        ("1002_AU-X-1", None),
        ("AU-X-1", None),  # one of theirs
        ("1001", "already exists in the cohort pedigree"),
        ("has space", "no space, comma or slash"),
        ("a,b", "no space, comma or slash"),
    ],
)
def test_the_target_id(data_dir, target, error):
    merge = merge_of(data_dir, "1002", "AU-X-1")
    found = [p.message for p in pm.target_problems(merge, target)]
    assert (found == []) if error is None else any(error in m for m in found)


def test_a_new_id_with_a_leftover_directory_is_refused(data_dir):
    merge = merge_of(data_dir, "1002", "AU-X-1")
    get_family_path(data_dir, "FAMX").mkdir(parents=True)
    (message,) = [p.message for p in pm.target_problems(merge, "FAMX")]
    assert "already exists: choose another ID" in message


def test_the_lines_are_merged_where_the_first_family_stood(data_dir):
    merge = merge_of(data_dir, "1002", "AU-X-1")
    draft = merge.original()
    draft.set_parent("S4", "father", "S5")  # S4 is S5's daughter after all
    text = demo(data_dir).read_text()
    new = pm.render_merge_text(text, demo(data_dir), merge, draft, "FAMX")
    rows = [
        line for line in new.splitlines() if line.startswith(("1002", "AU-X-1", "FAMX"))
    ]
    assert rows == [
        "FAMX\tS4\tS5\t0\t2.0\t-9.0",
        "FAMX\tS5\t-9\t\t1\t2",
        "FAMX\tS6\tS5\t0\t2\t1",
    ]
    # Into one of theirs: its unchanged lines are kept as they were.
    keep = pm.render_merge_text(text, demo(data_dir), merge, merge.original(), "AU-X-1")
    assert (
        "AU-X-1\tS5\t-9\t\t1\t2\n" in keep and "AU-X-1\tS4\t0\t0\t2.0\t-9.0\n" in keep
    )
    assert "1001\tS1" in keep


def test_a_merge_under_a_new_id(data_dir, curator):
    notes = data_dir / "notes" / "notes.tsv"
    append_tsv_row(notes, ["FID", "Sample", "Message", "User", "Timestamp"],
                   ["AU-X-1", "S6", "seen at the clinic", "alice", "2026-11-01"])  # fmt: skip
    acks_before = (data_dir / "validations" / "sampleqc.tsv").read_text()
    merge = merge_of(data_dir, "1002", "AU-X-1")
    draft = merge.original()
    records = pm.save_merge(
        merge, draft, "1002_AU-X-1", user="alice", comment="one family", now=NOW
    )
    assert [(r.family_id, r.rerun) for r in records] == [
        ("1002", pm.RERUN_MERGED),
        ("AU-X-1", pm.RERUN_MERGED),
        ("1002_AU-X-1", pe.RERUN_NORMAL),
    ]
    target = records[-1]
    assert target.changes[0] == "1002, AU-X-1 merged into 1002_AU-X-1 (a new family)"
    assert (
        "curation re-keyed: 1 in validations/sampleqc.tsv, 1 in notes/notes.tsv"
        in target.changes
    )
    assert records[0].changes == ("merged into 1002_AU-X-1",)

    store = DataStore(data_dir)
    store.load()
    families = store.cohorts["DEMO"].families
    assert "1002" not in families and "AU-X-1" not in families
    assert [s.sample_id for s in families["1002_AU-X-1"].samples] == ["S4", "S5", "S6"]
    family_file = get_family_path(data_dir, "1002_AU-X-1") / "1002_AU-X-1.pedigree.tsv"
    assert family_file.read_text().count("1002_AU-X-1\t") == 3

    # Curation re-keyed; the cross-family acknowledgement left as it was.
    acks = (data_dir / "validations" / "sampleqc.tsv").read_text()
    assert "SEX_MISMATCH\t1002_AU-X-1\tS6" in acks
    assert "CRYPTIC_RELATED\t1001,AU-X-1\tS2\tS5" in acks
    assert "1002_AU-X-1\tS6\tseen at the clinic" in notes.read_text()
    backup = (
        pe.edits_dir(data_dir, "1002_AU-X-1") / f"1002_AU-X-1.{target.edit_id}.curation"
    )
    assert (backup / "validations" / "sampleqc.tsv").read_text() == acks_before

    # The history: the target has every family's rows before; each old family
    # says where it went; the index has the three.
    before = pe.edit_rows(data_dir, "1002_AU-X-1", target.edit_id, "before")
    assert "1002\tS4" in before and "AU-X-1\tS6" in before
    assert pe.read_family_edits(data_dir, "AU-X-1")[-1].changes == (
        "merged into 1002_AU-X-1",
    )
    orphans = pm.orphaned_families(data_dir, pe.read_edit_index(data_dir), families)
    assert [(o.family_id, o.into) for o in orphans] == [
        ("1002", "1002_AU-X-1"), ("AU-X-1", "1002_AU-X-1"),
    ]  # fmt: skip


def test_a_merge_into_one_of_the_families(data_dir, curator):
    merge = merge_of(data_dir, "1002", "AU-X-1")
    records = pm.save_merge(
        merge,
        merge.original(),
        "AU-X-1",
        user="alice",
        comment="S4 belongs here",
        now=NOW,
    )
    assert [(r.family_id, r.rerun) for r in records] == [
        ("1002", pm.RERUN_MERGED),
        ("AU-X-1", pe.RERUN_CLEAN),
    ]
    assert records[-1].changes[0] == "1002 merged into AU-X-1"
    text = demo(data_dir).read_text()
    assert "1002\t" not in text and "AU-X-1\tS4\t0\t0\t2.0\t-9.0" in text
    # Its own acknowledgement is not re-keyed (it is the target already).
    assert "curation re-keyed" not in " ".join(records[-1].changes)


def test_a_refused_merge_changes_nothing(data_dir, curator):
    before = demo(data_dir).read_text()
    merge = merge_of(data_dir, "1002", "AU-X-1")
    with pytest.raises(pe.PedigreeEditError, match="already exists"):
        pm.save_merge(
            merge, merge.original(), "1001", user="alice", comment="x", now=NOW
        )
    with pytest.raises(pe.PedigreeEditError, match="comment is required"):
        pm.save_merge(
            merge, merge.original(), "FAMX", user="alice", comment=" ", now=NOW
        )
    assert demo(data_dir).read_text() == before
    assert not pe.index_file(data_dir).exists()


def test_a_viewer_cannot_merge(data_dir, monkeypatch):
    monkeypatch.setattr(pe, "can_write", lambda: False)
    merge = merge_of(data_dir, "1002", "AU-X-1")
    with pytest.raises(PermissionError):
        pm.save_merge(merge, merge.original(), "FAMX", user="eve", comment="x", now=NOW)
