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
