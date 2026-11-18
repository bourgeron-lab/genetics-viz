"""Family status: the record, the reasons to look again, the yield."""

import os
from pathlib import Path

import datagen
import pytest

from genetics_viz.utils import family_status as fs
from genetics_viz.utils import pedigree_edit as pe
from genetics_viz.utils import pedigree_merge as pm
from genetics_viz.utils.tsv import append_tsv_row


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    return datagen.build_data_dir(tmp_path)


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(fs, "can_write", lambda: True)


def _edit(timestamp: str, awaiting: bool, rerun: str = pe.RERUN_NORMAL):
    record = pe.EditRecord(
        "1001",
        "E1",
        timestamp,
        "bob",
        "father fixed",
        ("S3: father",),
        ("DEMO",),
        rerun,
    )
    return pe.EditStatus(record, awaiting)


def _record(status: str, timestamp: str = "2026-11-01T10:00:00") -> fs.StatusRecord:
    return fs.StatusRecord("1001", status, "", "alice", timestamp)


def test_every_status_is_worded_and_counted():
    assert list(fs.STATUSES)[0] == fs.NOT_REVIEWED
    for key, info in fs.STATUSES.items():
        assert info["label"] and info["color"].startswith("#"), key
        assert info["yield"] in ("solved", "reviewed", "open", "excluded"), key
    assert fs.SOLVED == {"solved", "partially_solved"}
    assert fs.EXCLUDED == {"not_analysable"}
    assert fs.status_info("bogus")["label"] == "bogus"


def test_save_last_row_wins_and_history(data_dir, curator):
    fs.save_family_status(
        data_dir, "1001", "in_review", "", user="alice", timestamp="2026-11-01T10:00:00"
    )
    fs.save_family_status(
        data_dir, "1001", "solved", "SHANK3\tde novo\nfits", user="bob",
        timestamp="2026-11-02T10:00:00",
    )  # fmt: skip
    record = fs.family_statuses(data_dir)["1001"]
    assert (record.status, record.user, record.date) == ("solved", "bob", "2026-11-02")
    assert record.comment == "SHANK3 de novo fits"
    assert [r.status for r in fs.family_history(data_dir, "1001")] == [
        "solved",
        "in_review",
    ]
    assert "1002" not in fs.family_statuses(data_dir)


def test_save_is_refused_to_viewers_and_unknown_statuses(data_dir, monkeypatch):
    monkeypatch.setattr(fs, "can_write", lambda: False)
    with pytest.raises(PermissionError):
        fs.save_family_status(data_dir, "1001", "solved", user="eve")
    monkeypatch.setattr(fs, "can_write", lambda: True)
    with pytest.raises(ValueError):
        fs.save_family_status(data_dir, "1001", "done", user="alice")
    assert not fs.status_file(data_dir).exists()


def test_no_record_is_not_reviewed_and_has_no_date_reason():
    state = fs.family_state("1001", None, edit=_edit("2026-11-10T10:00:00", True))
    assert state.status == fs.NOT_REVIEWED
    assert not state.reasons and not state.reopened


@pytest.mark.parametrize(
    "edited, awaiting, rerun, code",
    [
        ("2026-11-10T10:00:00", True, pe.RERUN_NORMAL, "pedigree_waiting"),
        ("2026-11-10T10:00:00", False, pe.RERUN_NORMAL, "pedigree_ready"),
        # A phenotype edit reruns nothing, but changes who is affected.
        ("2026-11-10T10:00:00", False, pe.RERUN_NONE, "pedigree_ready"),
        ("2026-10-10T10:00:00", False, pe.RERUN_NORMAL, None),
    ],
)
def test_a_later_pedigree_edit_reopens(edited, awaiting, rerun, code):
    state = fs.family_state(
        "1001",
        _record("unsolved"),
        edit=_edit(edited, awaiting, rerun),
    )
    assert [r.code for r in state.reasons] == ([code] if code else [])
    assert state.reopened is bool(code)
    assert state.waiting_rerun is (code == "pedigree_waiting")
    assert state.status == "unsolved"


def test_diagnostics_against_the_status():
    pathogenic = fs.FamilyDiagnostics(frozenset({"pathogenic"}), "2026-10-01T10:00:00")
    later = fs.FamilyDiagnostics(frozenset({"uncertain"}), "2026-11-05T10:00:00")

    def codes(record, diag):
        return [
            r.code for r in fs.family_state("1001", record, diagnostics=diag).reasons
        ]

    assert codes(_record("solved"), pathogenic) == []
    assert codes(_record("solved"), None) == ["solved_without_diagnostic"]
    assert codes(_record("unsolved"), pathogenic) == ["diagnostic_not_solved"]
    # A family nobody reviewed yet, with a pathogenic diagnostic: set a status.
    assert codes(None, pathogenic) == ["diagnostic_not_solved"]
    assert codes(_record("in_review"), pathogenic) == []
    assert codes(_record("not_analysable"), pathogenic) == []
    assert codes(_record("candidate"), later) == ["diagnostic_after"]


def test_qc_and_results_written_after_the_status():
    state = fs.family_state(
        "1001",
        _record("unsolved"),
        qc_high_written="2026-11-03T10:00:00",
        results=[("sampleqc", "2026-11-04T10:00:00"), ("wombat", "2026-10-01")],
    )
    assert [r.code for r in state.reasons] == ["qc_after", "results_after"]
    assert "sampleqc" in state.reasons[1].text and "wombat" not in state.reasons[1].text
    assert not state.reopened
    # Earlier writes, or a st_mtime before the status, say nothing.
    old = fs.as_moment("2026-10-01T00:00:00").timestamp()
    assert not fs.family_state(
        "1001", _record("unsolved"), qc_high_written=old, results=[("upd", old)]
    ).reasons


def test_yield_counts_analysable_families():
    def st(fid, status, affected=True, **kw):
        record = fs.StatusRecord(fid, status, "", "a", "2026-11-01") if status else None
        return fs.family_state(fid, record, affected=affected, **kw)

    states = [
        st("A", "solved", diagnostics=fs.FamilyDiagnostics(frozenset({"pathogenic"}))),
        st("B", "partially_solved"),
        st("C", "unsolved"),
        st("D", "candidate"),
        st("E", "in_review"),
        st("F", None),
        st("G", "not_analysable"),
        # No affected member and no status: not in the yield.
        st("H", None, affected=False),
        st("I", "unsolved", edit=_edit("2026-11-10T10:00:00", True)),
    ]
    counts = fs.count_statuses(states)
    assert counts.total == 9 and counts.excluded == 2
    assert counts.solved == 2 and counts.reviewed == 5 and counts.analysable == 7
    assert counts.not_reviewed == 1
    assert counts.reopened == 1
    assert counts.to_check == 2  # B (no pathogenic) and I (reopened)
    assert counts.yield_text() == "Solved 2 of 5 reviewed (40%) · 7 analysable"
    assert fs.count_statuses([]).yield_text() == "No family to review"


def test_badges_and_filter_categories():
    state = fs.family_state(
        "1001",
        fs.StatusRecord("1001", "unsolved", "non-paternity", "alice", "2026-11-01"),
        edit=_edit("2026-11-10T10:00:00", True),
    )
    badges = fs.status_badges(state, curator=False)
    assert [b["label"] for b in badges] == ["Unsolved", "Reopened"]
    assert "non-paternity" not in badges[0]["tooltip"]
    assert "non-paternity" in fs.status_badges(state, curator=True)[0]["tooltip"]
    assert fs.filter_categories(state) == ["Unsolved", fs.REOPENED, fs.TO_CHECK]
    assert set(fs.filter_categories(state)) <= set(fs.filter_options())
    lonely = fs.family_state("2001", None, affected=False)
    assert fs.status_badges(lonely, False)[0]["label"] == "Not reviewed (no affected)"


def test_diagnostics_by_family_skip_ignored_rows(data_dir):
    header = datagen.DIAGNOSTIC_COLUMNS
    path = data_dir / "diagnostics" / "snvs.tsv"
    for row in datagen.REPORT_DIAGNOSTICS:
        append_tsv_row(data_dir / "diagnostics" / row[0], header, row[1])
    diag = fs.diagnostics_by_family(data_dir)["1001"]
    assert diag.verdicts == {"pathogenic", "uncertain", "benign"}
    # The ignored pathogenic row of 2026-10-04 is not the latest.
    assert diag.latest == "2026-10-03T10:00:00"
    assert diag.solving
    assert path.exists()


def test_load_family_state_reads_everything(data_dir, curator):
    store_family = _family(data_dir, "1001")
    fs.save_family_status(
        data_dir, "1001", "unsolved", user="alice", timestamp="2000-01-01T00:00:00"
    )
    append_tsv_row(
        pe.index_file(data_dir),
        pe.INDEX_HEADER,
        ["1001", "E1", "2026-11-10T10:00:00", "bob", "fix", "S3: father", "DEMO",
         pe.RERUN_NORMAL],
    )  # fmt: skip
    state = fs.load_family_state(data_dir, "1001", store_family, [_cohort(data_dir)])
    codes = [r.code for r in state.reasons]
    # DEMO's last run started on 2026-10-01: the edit waits for one.
    assert codes[0] == "pedigree_waiting"
    # 1001 has an open NOT_PARENT (high) in its QC tables, written today.
    assert "qc_after" in codes
    assert "results_after" in codes
    assert state.affected


def test_the_status_file_is_re_keyed_by_a_merge(data_dir, curator, monkeypatch):
    monkeypatch.setattr(pe, "can_write", lambda: True)
    assert fs.status_file(data_dir) in pm.curation_files(data_dir)
    fs.save_family_status(
        data_dir, "AU-X-1", "unsolved", user="alice", timestamp="2026-11-01T10:00:00"
    )
    demo = data_dir / "cohorts" / "DEMO" / "DEMO.pedigree.tsv"
    merge = pm.load_merge(data_dir, ["1002", "AU-X-1"], [[demo], [demo]])
    pm.save_merge(merge, merge.original(), "1002_AU-X-1", user="alice", comment="one")
    statuses = fs.family_statuses(data_dir)
    assert "AU-X-1" not in statuses
    assert statuses["1002_AU-X-1"].status == "unsolved"
    # The merge is an edit after the status: the family opens reopened.
    edits = fs.cohort_edits(None, data_dir)
    state = fs.family_state(
        "1002_AU-X-1", statuses["1002_AU-X-1"], edit=edits["1002_AU-X-1"]
    )
    assert state.reopened


def test_newest_and_step_outputs(data_dir):
    outputs = fs.step_outputs(data_dir, "1001")
    assert outputs and all(isinstance(m, float) for _, m in outputs)
    path = data_dir / "nothing.tsv"
    assert fs.newest([path]) is None
    path.write_text("x")
    os.utime(path, (1_000_000, 1_000_000))
    assert fs.newest([path, data_dir / "nope"]) == 1_000_000


def _cohort(data_dir: Path):
    from genetics_viz.models import DataStore

    store = DataStore(data_dir)
    store.load()
    return store.cohorts["DEMO"]


def _family(data_dir: Path, fid: str):
    return _cohort(data_dir).families[fid]
