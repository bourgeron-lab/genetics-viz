"""Run requests: the files genetics-viz writes and swiftbar answers (utils/run_requests)."""

import json
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace

import pytest

from genetics_viz.utils import cohort_state as cs
from genetics_viz.utils import run_requests as rr


@pytest.fixture
def cohort(tmp_path, monkeypatch):
    """A cohort directory whose params ask for four steps; the caller is a curator."""
    monkeypatch.setattr(rr, "can_write", lambda: True)
    monkeypatch.setattr(rr, "is_admin", lambda: False)
    d = tmp_path / "cohorts" / "C"
    d.mkdir(parents=True)
    (d / "C.pedigree.tsv").write_text("F1\tS1\t0\t0\t1\t2\n")
    params = d / "C.params.yml"
    params.write_text("steps: [sampleqc, genomeqc, 'wombat,upd']\n")
    return SimpleNamespace(
        name="C", pedigree_file=d / "C.pedigree.tsv", params_file=params
    )


def folder(cohort):
    return rr.requests_dir(cohort.pedigree_file.parent)


def handle(cohort, rid, **content):
    (folder(cohort) / f"{rid}.handled.json").write_text(json.dumps(content))


def test_allowed_steps_split_a_comma_string(cohort):
    assert rr.allowed_steps(cohort) == ["sampleqc", "genomeqc", "wombat", "upd"]
    assert rr.allowed_steps(SimpleNamespace(params_file=None)) == []


def test_a_request_is_written_whole_and_read_back(cohort):
    req = rr.write_request(cohort, ["sampleqc", "genomeqc"], reason=" edits ", user="ana",
                           families=["1001"], origin="family 1001, Sample QC tab")  # fmt: skip
    assert rr.ID_PATTERN.match(req.id) and req.state == rr.OPEN
    raw = json.loads((folder(cohort) / f"{req.id}.json").read_text())
    assert raw["schema_version"] == 1 and raw["cohort"] == "C"
    assert (
        raw["steps"] == ["sampleqc", "genomeqc"]
        and raw["clean_stale_families"] is False
    )
    assert raw["reason"] == "edits" and raw["requested_by"] == "ana"
    assert raw["families"] == ["1001"]
    # Nothing but the request and the lock: no temporary file left behind
    assert sorted(p.name for p in folder(cohort).iterdir()) == [
        ".lock",
        f"{req.id}.json",
    ]
    assert (
        rr.command_line("C", req.steps, False)
        == "ghfc-ngs C --submit --steps sampleqc,genomeqc"
    )


def test_a_normal_run_and_a_clean(cohort):
    req = rr.write_request(
        cohort, None, clean=True, reason="members changed", user="ana"
    )
    assert req.steps is None and req.clean_stale_families
    assert rr.steps_text(req.steps) == "normal run"
    assert (
        rr.command_line("C", None, True) == "ghfc-ngs C --submit --clean-stale-families"
    )


def test_refusals(cohort, monkeypatch):
    with pytest.raises(ValueError, match="why"):
        rr.write_request(cohort, None, reason="  ", user="ana")
    with pytest.raises(ValueError, match="Not a step of C's params file: roh"):
        rr.write_request(cohort, ["sampleqc", "roh"], reason="x", user="ana")
    with pytest.raises(ValueError, match="at least one step"):
        rr.write_request(cohort, [], reason="x", user="ana")
    monkeypatch.setattr(rr, "can_write", lambda: False)
    with pytest.raises(PermissionError):
        rr.write_request(cohort, None, reason="x", user="ana")
    assert not folder(cohort).exists() or not list(folder(cohort).glob("*.json"))


def test_one_open_request_per_user_and_cohort(cohort):
    first = rr.write_request(cohort, ["sampleqc"], reason="a", user="ana")
    with pytest.raises(ValueError, match="Withdraw it first"):
        rr.write_request(cohort, ["wombat"], reason="b", user="ana")
    other = rr.write_request(
        cohort, ["wombat"], reason="b", user="ben"
    )  # another curator
    rr.withdraw_request(cohort.pedigree_file.parent, first.id, user="ana")
    again = rr.write_request(cohort, ["wombat"], reason="c", user="ana")
    states = {r.id: r.state for r in rr.read_requests(cohort.pedigree_file.parent)}
    assert states == {first.id: rr.WITHDRAWN, other.id: rr.OPEN, again.id: rr.OPEN}
    mine = rr.open_request_of(rr.read_requests(cohort.pedigree_file.parent), "ana")
    assert mine.id == again.id


def test_withdraw_is_for_the_author_or_an_administrator(cohort, monkeypatch):
    req = rr.write_request(cohort, None, reason="a", user="ana")
    d = cohort.pedigree_file.parent
    monkeypatch.setattr(rr, "is_admin", lambda: False)
    with pytest.raises(PermissionError):
        rr.withdraw_request(d, req.id, user="ben")
    monkeypatch.setattr(rr, "is_admin", lambda: True)
    rr.withdraw_request(d, req.id, user="admin")
    (read,) = rr.read_requests(d)
    assert read.state == rr.WITHDRAWN and read.withdrawn_by == "admin"
    with pytest.raises(ValueError, match="already withdrawn"):
        rr.withdraw_request(d, req.id, user="admin")
    with pytest.raises(ValueError, match="Not a request"):
        rr.withdraw_request(d, "../../etc/passwd", user="admin")


def test_handled_requests_and_their_run(cohort):
    d = cohort.pedigree_file.parent
    started = rr.write_request(cohort, ["sampleqc", "genomeqc"], reason="a", user="ana")
    dismissed = rr.write_request(cohort, None, reason="b", user="ben")
    handle(cohort, started.id, status="started", user="fred", at="2026-11-25T09:00:00+01:00",
           job_id="48945575", steps_sent=["sampleqc"], clean_sent=False)  # fmt: skip
    handle(cohort, dismissed.id, status="dismissed", user="fred", note="ran yesterday")
    by_id = {r.id: r for r in rr.read_requests(d)}
    assert by_id[started.id].state == rr.STARTED and by_id[started.id].edited
    assert by_id[started.id].handled.steps_sent == ("sampleqc",)
    assert by_id[started.id].outcome is None  # no record names the job yet
    assert by_id[dismissed.id].state == rr.DISMISSED
    assert by_id[dismissed.id].handled.note == "ran yesterday"
    # The run record names the job: its status is the request's outcome
    record = {
        "status": "success",
        "run_id": "r",
        "launch": {"slurm_job_id": "48945575"},
    }
    (d / cs.STATE_FILENAME).write_text(json.dumps(
        {"schema_version": 1, "last_run": record, "history": [record]}))  # fmt: skip
    cs.clear_state_cache()
    assert {r.id: r.outcome for r in rr.read_requests(d)}[started.id] == "success"


def test_files_that_are_not_requests_are_skipped(cohort):
    f = folder(cohort)
    f.mkdir()
    (f / "20261124T101500-abcdef.json").write_text("{not json")
    (f / "nope.json").write_text(json.dumps({"schema_version": 1}))
    (f / "20261124T101501-abcdef.json").write_text(json.dumps({"schema_version": 2}))
    (f / "20261124T101502-abcdef.json").write_text(json.dumps(
        {"schema_version": 1, "steps": ["a", 3, None], "reason": "x" * 2000,
         "families": list(range(100)), "clean_stale_families": "yes"}))  # fmt: skip
    (req,) = rr.read_requests(cohort.pedigree_file.parent)
    assert req.steps == ("a",) and len(req.reason) == rr.REASON_MAX
    assert len(req.families) == rr.FAMILIES_MAX and req.clean_stale_families is False


def test_recent_keeps_open_and_recently_settled(cohort):
    d = cohort.pedigree_file.parent
    old = rr.write_request(cohort, None, reason="a", user="ana")
    handle(cohort, old.id, status="dismissed", at="2026-01-01T00:00:00+00:00")
    fresh = rr.write_request(cohort, None, reason="b", user="ana")
    now = datetime(2026, 11, 25, tzinfo=timezone.utc)
    assert [r.id for r in rr.recent(rr.read_requests(d), now)] == [fresh.id]
    later = now + timedelta(days=400)
    assert [r.id for r in rr.recent(rr.read_requests(d), later)] == [
        fresh.id
    ]  # still open


def test_requests_are_off_unless_the_directory_turns_them_on(tmp_path):
    assert not rr.requests_enabled(tmp_path)  # no config loaded


def test_a_familys_pending_steps_become_the_request(cohort):
    from genetics_viz.components.pedigree_editor import _request_steps

    # Wombat brings its cohort merge when the params run it
    cohort.params_file.write_text(
        "steps: [sampleqc, genomeqc, wombat, wombat_cohort, upd]\n"
    )
    assert _request_steps(cohort, {"wombat", "sampleqc"}) == [
        "sampleqc", "wombat", "wombat_cohort"
    ]  # fmt: skip
    # None of the pending steps in the params: a normal run
    assert _request_steps(cohort, {"expansionhunter"}) is None
