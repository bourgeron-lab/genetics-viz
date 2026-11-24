"""Pipeline run record, parameters file and the fetched reference (cohort_state, pipeline_params, ghfc_ngs_docs)."""

import json
import urllib.error
from datetime import datetime, timedelta, timezone

import pytest

from genetics_viz.utils import cohort_state as cs
from genetics_viz.utils import ghfc_ngs_docs as docs
from genetics_viz.utils import pipeline_params as pp


def write_state(tmp_path, payload):
    (tmp_path / cs.STATE_FILENAME).write_text(
        payload if isinstance(payload, str) else json.dumps(payload)
    )
    cs.clear_state_cache()
    return cs.read_cohort_state(tmp_path)


def test_synthetic_record_parses(data_dir):
    state = cs.read_cohort_state(data_dir / "cohorts" / "DEMO")
    run = state.last_run
    assert not state.unsupported and run.status == "success"
    assert (
        run.duration == timedelta(hours=2, minutes=30)
        and run.short_commit == "abcdef12"
    )
    assert run.completion["ancestry"].complete and not run.completion["roh"].measured


def test_absent_or_broken_records(tmp_path):
    assert cs.read_cohort_state(tmp_path) is None
    assert write_state(tmp_path, "{not json") is None
    assert write_state(tmp_path, "[1, 2]") is None


@pytest.mark.parametrize("version", [2, "1", True, None])
def test_unknown_schema_is_refused_whole(tmp_path, version):
    state = write_state(
        tmp_path, {"schema_version": version, "last_run": {"status": "success"}}
    )
    assert state.unsupported and state.last_run is None
    assert cs.status_meta(state)["label"] == "unreadable record"


def test_null_completion_and_bool_counts(tmp_path):
    state = write_state(
        tmp_path,
        {
            "schema_version": 1,
            "last_run": {
                "status": "failed",
                "completion": {"a": {"done": True, "total": 3}, "b": None},
            },
        },
    )
    completion = state.last_run.completion
    assert completion["a"].done is None and not completion["b"].measured


def test_stale_running_record():
    old = datetime.now(timezone.utc) - timedelta(days=2)
    assert cs.RunRecord(status="running", started_at=old).is_stale_running
    assert cs.RunRecord(
        status="running", started_at=old.replace(tzinfo=None)
    ).is_stale_running
    assert not cs.RunRecord(
        status="running", started_at=datetime.now(timezone.utc)
    ).is_stale_running
    assert (
        cs.record_status_meta(cs.RunRecord(status="running", started_at=old))["label"]
        == "running (stale)"
    )
    assert (
        cs.record_status_meta(cs.RunRecord(status="brand_new"))["label"] == "brand_new"
    )


def test_file_drift(tmp_path):
    f = tmp_path / "p.tsv"
    f.write_text("a")
    sha = cs.sha256_file(f)
    assert cs.file_drift(sha, f) == "unchanged"
    assert cs.file_drift("0" * 64, f) == "changed"
    assert cs.file_drift(None, f) == cs.file_drift(sha, tmp_path / "gone") == "unknown"


def test_step_rows_without_a_record_follow_the_params():
    rows = cs.build_step_rows(None, ["roh", "ancestry", "zzz_new_step"])
    assert [r.step for r in rows][-1] == "zzz_new_step"  # unknown steps appended
    assert all(not r.measured and r.in_params for r in rows)
    assert cs.summarize_steps(rows).text == "3 steps, none measured"
    assert rows[0].detail() == "not measured"


def test_step_rows_with_a_record(data_dir):
    state = cs.read_cohort_state(data_dir / "cohorts" / "DEMO")
    rows = {r.step: r for r in cs.build_step_rows(state, ["ancestry", "roh"])}
    assert rows["ancestry"].requested and rows["ancestry"].detail() == "3/3 families"
    assert rows["ancestry"].detail(compact=True) == "3/3"
    summary = cs.summarize_steps(list(rows.values()))
    assert summary.text == "1/1 steps complete, 1 not measured" and summary.all_complete


@pytest.mark.parametrize(
    "seconds, text",
    [(59, "59s"), (61, "1m 01s"), (3 * 3600, "3h 00m"), (-5, "—"), (None, "—")],
)
def test_format_duration(seconds, text):
    assert (
        cs.format_duration(None if seconds is None else timedelta(seconds=seconds))
        == text
    )


# --- params file ----------------------------------------------------------------


def test_params_file_discovery_and_parsing(tmp_path):
    assert pp.find_params_file(tmp_path, "C") is None
    (tmp_path / "C.params.yaml").write_text("steps: ancestry\n")
    (tmp_path / "C.params.yml").write_text(
        "# comment kept\nsteps: [ancestry, null, 3]\n"
    )
    path = pp.find_params_file(tmp_path, "C")
    assert path.name == "C.params.yml"  # .yml preferred
    assert pp.get_requested_steps(pp.load_params(path)) == ["ancestry", "3"]
    assert pp.read_params_text(path).startswith("# comment kept")
    assert pp.get_requested_steps({"steps": "roh"}) == ["roh"]
    assert pp.get_requested_steps({}) == []


def test_malformed_params_are_empty(tmp_path):
    (tmp_path / "bad.yml").write_text("steps: [unclosed\n")
    (tmp_path / "list.yml").write_text("- a\n- b\n")
    assert pp.load_params(tmp_path / "bad.yml") == {}
    assert pp.load_params(tmp_path / "list.yml") == {}


# --- ghfc-ngs reference ------------------------------------------------------------


@pytest.mark.parametrize(
    "markdown, expected",
    [
        ("[top](#steps)", "[top](#steps)"),
        ("[web](https://x.org)", "[web](https://x.org)"),
        (
            "[r](../README.md#install)",
            f'<a href="{docs._BLOB}README.md#install" target="_blank" rel="noopener">r</a>',
        ),
        (
            "[ex](../params_example/)",
            f'<a href="{docs._TREE}params_example/" target="_blank" rel="noopener">ex</a>',
        ),
        (
            "[`cfg`](../nextflow.config)",
            f'<a href="{docs._BLOB}nextflow.config" target="_blank" rel="noopener"><code>cfg</code></a>',
        ),
    ],
)
def test_absolutise_links(markdown, expected):
    assert docs.absolutise_links(markdown) == expected


class FakeResponse:
    def __init__(self, body: bytes):
        self.body = body

    def read(self):
        return self.body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_fetch_is_cached_and_failures_are_reported(monkeypatch):
    calls = []

    def ok(request, timeout):
        calls.append(request.full_url)
        return FakeResponse(b"# Params\n[r](../README.md)")

    monkeypatch.setattr(docs.urllib.request, "urlopen", ok)
    first = docs.fetch_params_doc()
    assert first.ok and "github.com" in first.markdown
    assert docs.fetch_params_doc() is first and len(calls) == 1
    docs.fetch_params_doc(force=True)
    assert len(calls) == 2

    def fail(request, timeout):
        raise urllib.error.URLError("offline")

    monkeypatch.setattr(docs.urllib.request, "urlopen", fail)
    failed = docs.fetch_params_doc(force=True)
    assert not failed.ok and "offline" in failed.error


def test_held_families(data_dir):
    import datagen

    state = cs.read_cohort_state(data_dir / "cohorts" / "DEMO")
    assert state.last_run.held_families is None  # not recorded
    datagen.add_held_families(data_dir)
    cs.clear_state_cache()
    (held,) = cs.read_cohort_state(data_dir / "cohorts" / "DEMO").last_run.held_families
    assert (held.family, held.cohorts) == ("1001", ("NOHEAD",))
    assert held.describe() == [
        "differs from NOHEAD: S3: father S1 -> 0, mother S2 -> 0"
    ]


def test_held_families_parse_defensively():
    assert cs._parse_held(None) is None
    assert cs._parse_held([]) == []
    held = cs._parse_held(
        [{"family": "F"}, {"cohorts": ["A"]}, "junk", {"family": "G", "cohorts": ["A", 3],
          "changes": {"A": ["x", 1], "B": "y"}}]
    )  # fmt: skip
    assert [(h.family, h.cohorts, h.changes) for h in held] == [
        ("F", (), {}),
        ("G", ("A",), {"A": ("x",)}),
    ]
    assert held[0].describe() == []


# --- per-step record (last_run_by_step) ----------------------------------------------


def _record(run_id, status, started, steps, *, ped="p1", par="q1", done=None):
    """A run record asked for ``steps``, each measured ``done``/10 (``None``: not measured)."""
    every = [
        "alignment",
        "wombat",
        "sampleqc",
        "genomeqc",
        "upd",
        "denovo",
        "expansionhunter",
    ]
    return {
        "run_id": run_id, "status": status, "started_at": started,
        "finished_at": None if status == "running" else started,
        "steps_requested": steps,
        "pedigree": {"sha256": ped}, "params_file": {"sha256": par},
        "completion_measured": "before" if status == "running" else "after",
        "completion": {s: ({"done": done if done is not None else 10, "total": 10,
                            "pct": (done if done is not None else 10) * 10.0} if s in steps else None)
                       for s in every},
    }  # fmt: skip


FULL = [
    "alignment",
    "wombat",
    "sampleqc",
    "genomeqc",
    "upd",
    "denovo",
    "expansionhunter",
]
QC = ["sampleqc", "genomeqc"]


def _entry(record, step):
    keys = ("run_id", "status", "started_at", "finished_at", "completion_measured")
    return {**{k: record[k] for k in keys}, "outputs_may_be_incomplete": False,
            "pedigree_sha256": record["pedigree"]["sha256"],
            "params_sha256": record["params_file"]["sha256"],
            "completion": record["completion"][step]}  # fmt: skip


def _file(full, qc, *, recorded):
    payload = {"schema_version": 1, "last_run": qc, "last_successful_run": qc,
               "history": [qc, full]}  # fmt: skip
    if recorded:
        payload["last_run_by_step"] = {
            s: {"last": _entry(qc if s in QC else full, s),
                "last_success": _entry(qc if s in QC else full, s)}
            for s in FULL
        }  # fmt: skip
    return payload


@pytest.mark.parametrize("recorded", [True, False])
def test_a_qc_only_run_leaves_the_other_steps_their_own_run(tmp_path, recorded):
    full = _record("r1", "success", "2026-11-19T10:00:00+01:00", FULL, done=9)
    qc = _record("r2", "success", "2026-11-25T09:00:00+01:00", QC, ped="p2")
    state = write_state(tmp_path, _file(full, qc, recorded=recorded))
    assert state.by_step_recorded is recorded
    rows = {r.step: r for r in cs.build_step_rows(state, FULL)}
    assert rows["sampleqc"].latest and rows["sampleqc"].detail() == "10/10 samples"
    wombat = rows["wombat"]
    assert wombat.requested and not wombat.latest and wombat.completion.done == 9
    assert wombat.run_at == datetime.fromisoformat("2026-11-19T10:00:00+01:00")
    assert cs.summarize_steps(list(rows.values())).text == "2/7 steps complete"


def test_derived_record_follows_the_pipeline_rules(tmp_path):
    ok = _record("r1", "success", "2026-11-01T10:00:00+00:00", ["wombat", "upd"])
    failed = _record("r2", "failed", "2026-11-02T10:00:00+00:00", ["wombat", "bogus"])
    invalid = {
        **_record("r3", "failed", "2026-11-03T10:00:00+00:00", ["upd"]),
        "completion": None,
    }
    state = write_state(tmp_path, {"schema_version": 1, "last_run": invalid,
                                   "last_successful_run": ok, "history": [invalid, failed]})  # fmt: skip
    # The last success (not in history) first, then the history oldest first;
    # the validation failure planned nothing, and "bogus" is no step.
    assert set(state.by_step) == {"wombat", "upd"}
    assert state.by_step["wombat"].last.run_id == "r2"
    assert state.by_step["wombat"].last_success.run_id == "r1"
    assert state.by_step["upd"].last.run_id == "r1"


def test_a_live_run_is_newer_than_any_finished_one(tmp_path):
    full = _record("r1", "success", "2026-11-19T10:00:00+00:00", FULL)
    now = datetime.now(timezone.utc).replace(microsecond=0)
    live = _record("r2", "running", now.isoformat(), QC, done=4)
    state = write_state(tmp_path, {"schema_version": 1, "last_run": live,
                                   "last_successful_run": full, "history": [full]})  # fmt: skip
    assert state.step_started("sampleqc") == now
    assert state.step_started("wombat") == datetime.fromisoformat(
        "2026-11-19T10:00:00+00:00"
    )
    rows = {r.step: r for r in cs.build_step_rows(state, FULL)}
    assert rows["sampleqc"].latest and rows["sampleqc"].completion.done == 4
    assert not rows["wombat"].latest


def test_a_pedigree_edit_waits_for_every_pedigree_step(tmp_path):
    full = _record("r1", "success", "2026-11-19T10:00:00+00:00", FULL)
    qc = _record("r2", "success", "2026-11-25T09:00:00+00:00", QC)
    state = write_state(tmp_path, _file(full, qc, recorded=True))
    # Wombat, UPD, de novo and STR last ran on the 19th: an edit of the 20th waits
    assert state.pedigree_steps_started(FULL) == datetime.fromisoformat(
        "2026-11-19T10:00:00+00:00"
    )
    # A cohort that runs only the QC steps has had them all since
    assert state.pedigree_steps_started(["alignment", *QC]) == datetime.fromisoformat(
        "2026-11-25T09:00:00+00:00"
    )
    # A pedigree step that never ran: every edit waits
    assert (
        state.pedigree_steps_started([*FULL, "roh"]) is not None
    )  # roh reads no snapshot
    assert cs.CohortState(path=tmp_path, by_step={}).pedigree_steps_started(QC) is None
    # No pedigree step in play: the last run's start, as before
    assert state.pedigree_steps_started(["alignment"]) == datetime.fromisoformat(
        "2026-11-25T09:00:00+00:00"
    )


def test_pedigree_run_started_reads_the_params_steps(tmp_path):
    from types import SimpleNamespace

    full = _record("r1", "success", "2026-11-19T10:00:00+00:00", FULL)
    qc = _record("r2", "success", "2026-11-25T09:00:00+00:00", QC)
    write_state(tmp_path, _file(full, qc, recorded=True))
    params = tmp_path / "C.params.yml"
    params.write_text("steps: [sampleqc, genomeqc]\n")
    cohort = SimpleNamespace(
        pedigree_file=tmp_path / "C.pedigree.tsv", params_file=params
    )
    assert cs.pedigree_run_started(cohort) == datetime.fromisoformat(
        "2026-11-25T09:00:00+00:00"
    )
    params.write_text("steps: [sampleqc, genomeqc, wombat]\n")
    assert cs.pedigree_run_started(cohort) == datetime.fromisoformat(
        "2026-11-19T10:00:00+00:00"
    )
    assert (
        cs.pedigree_run_started(SimpleNamespace(pedigree_file=tmp_path / "x" / "p"))
        is None
    )


def test_drift_is_per_step(tmp_path):
    pedigree = tmp_path / "C.pedigree.tsv"
    pedigree.write_text("edited\n")
    new = cs.sha256_file(pedigree)
    full = _record("r1", "success", "2026-11-19T10:00:00+00:00", FULL, ped="old")
    qc = _record("r2", "success", "2026-11-25T09:00:00+00:00", QC, ped=new)
    state = write_state(tmp_path, _file(full, qc, recorded=True))
    drift = cs.step_drift(state, pedigree)
    assert drift["sampleqc"] == drift["genomeqc"] == "unchanged"
    assert {s for s, d in drift.items() if d == "changed"} == set(FULL) - set(QC)
    assert set(cs.step_drift(state, None).values()) == {"unknown"}


def test_by_step_parses_defensively():
    assert cs._parse_by_step(None) is None and cs._parse_by_step([]) is None
    parsed = cs._parse_by_step(
        {
            "a": "junk",
            "b": {"last": None},
            "c": {"last": {"status": "failed", "completion": None}},
        }
    )
    assert list(parsed) == ["c"] and not parsed["c"].last.completion.measured
    assert parsed["c"].last_success is None
