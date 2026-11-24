"""Run requests on the cohort page: the Workflow list, the dialog, the shortcuts, who sees them."""

import asyncio
import json
from datetime import datetime
from pathlib import Path

import datagen
import pytest
import yaml
from ui_helpers import login, open_tab

REQUESTS = ".ghfc-ngs.requests"


@pytest.fixture
def config_file(config_file: Path) -> Path:
    """Run requests on for the synthetic data directory."""
    config = yaml.safe_load(config_file.read_text())
    config["data_directories"][0]["run_requests"] = True
    config_file.write_text(yaml.safe_dump(config))
    return config_file


def _stub_docs() -> None:
    import sys

    docs = sys.modules["genetics_viz.utils.ghfc_ngs_docs"]
    docs._cached = docs.FetchedDoc(
        markdown="# Stub reference", error=None, fetched_at=datetime.now()
    )


async def _workflow(user, username: str) -> None:
    await login(user, username)
    _stub_docs()
    await user.open("/cohort/DEMO?tab=workflow")
    await user.should_see("Pipeline status", retries=30)


async def wait_for(condition, tries: int = 100) -> None:
    for _ in range(tries):
        if condition():
            return
        await asyncio.sleep(0.05)
    raise AssertionError("condition never met")


def _requests(data_dir: Path) -> list:
    folder = data_dir / "cohorts" / "DEMO" / REQUESTS
    return sorted(p.name for p in folder.glob("*.json")) if folder.exists() else []


async def test_the_workflow_tab_lists_requests_and_what_became_of_them(user, data_dir):
    datagen.add_run_requests(data_dir)
    await _workflow(user, "curator")
    await user.should_see(marker="workflow-requests", retries=30)
    for rid in datagen.RUN_REQUESTS:
        await user.should_see(marker=f"run-request-{rid}", retries=30)
    await user.should_see("ROH after the new panel")
    # Started with fewer steps than asked, by the job the run record names
    await user.should_see(f"Slurm job {datagen.RUN_REQUEST_JOB}")
    await user.should_see("sent: Runs of homozygosity")
    await user.should_see("succeeded")
    await user.should_see("Nothing out of date")  # the dismissal's note
    # The admin's open request is not the curator's to withdraw
    await user.should_not_see(marker="run-request-list-withdraw-20261120T090000-open01")


async def test_a_curator_requests_a_run_then_withdraws_it(user, data_dir):
    await _workflow(user, "curator")
    await user.should_see(marker="run-request", retries=30)
    user.find(marker="run-request").click()
    await user.should_see(marker="run-request-dialog", retries=20)
    await user.should_see("ghfc-ngs DEMO --submit")  # a normal run by default
    user.find(marker="run-request-reason").type("Pedigree edits")
    user.find(marker="run-request-send").click()
    await wait_for(lambda: len(_requests(data_dir)) == 1)
    (name,) = _requests(data_dir)
    raw = json.loads((data_dir / "cohorts" / "DEMO" / REQUESTS / name).read_text())
    assert raw["requested_by"] == "curator" and raw["steps"] is None
    assert raw["reason"] == "Pedigree edits" and raw["origin"] == "Workflow tab"
    rid = name[: -len(".json")]
    # One open request per user: the button becomes the request, with Withdraw
    await user.should_see(marker="run-request-open", retries=30)
    await user.should_see(marker=f"run-request-{rid}", retries=30)
    await user.should_see(marker=f"run-request-list-withdraw-{rid}")
    user.find(marker=f"run-request-withdraw-{rid}").click()
    await wait_for(
        lambda: (
            data_dir / "cohorts" / "DEMO" / REQUESTS / f"{rid}.withdrawn.json"
        ).exists()
    )
    await user.should_see(marker="run-request", retries=30)


async def test_the_dialog_refuses_a_request_without_a_reason(user, data_dir):
    await _workflow(user, "curator")
    await user.should_see(marker="run-request", retries=30)
    user.find(marker="run-request").click()
    await user.should_see(marker="run-request-dialog", retries=20)
    user.find(marker="run-request-send").click()
    await user.should_see("Say why the run is needed", retries=20)
    assert _requests(data_dir) == []


async def test_the_workflow_button_counts_open_requests_for_curators(user, data_dir):
    datagen.add_run_requests(data_dir)
    await login(user, "curator")
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=20)

    def texts():
        try:
            return sorted(
                b.text for b in user.find(marker="workflow-tab-pills").elements
            )
        except AssertionError:
            return []

    await wait_for(lambda: "1" in texts())


@pytest.mark.parametrize("username", ["reader", "student"])
async def test_viewers_see_no_request(user, data_dir, username):
    datagen.add_run_requests(data_dir)
    await _workflow(user, username)
    await user.should_not_see(marker="workflow-requests")
    await user.should_not_see(marker="run-request")
    await user.should_not_see("ROH after the new panel")


async def test_the_cohort_sample_qc_tab_requests_the_qc_steps(user, data_dir):
    params = data_dir / "cohorts" / "DEMO" / "DEMO.params.yml"
    params.write_text(
        params.read_text().replace("  - roh\n", "  - roh\n  - sampleqc\n")
    )
    await login(user, "curator")
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=20)
    open_tab(user, "cohort-tabs", "Sample QC")
    await user.should_see(marker="sampleqc-request-run", retries=40)
    user.find(marker="sampleqc-request-run").click()
    await user.should_see(marker="run-request-dialog", retries=20)
    await user.should_see("ghfc-ngs DEMO --submit --steps sampleqc")
    user.find(marker="run-request-reason").type("New QC thresholds")
    user.find(marker="run-request-send").click()
    await wait_for(lambda: len(_requests(data_dir)) == 1)
    (name,) = _requests(data_dir)
    raw = json.loads((data_dir / "cohorts" / "DEMO" / REQUESTS / name).read_text())
    assert raw["steps"] == ["sampleqc"] and raw["origin"] == "cohort Sample QC tab"


async def test_an_administrator_turns_requests_off(user, config_file):
    await login(user, "admin")
    await user.open("/admin/directories")
    await user.should_see(marker="admin-run-requests-0", retries=20)
    switch = user.find(marker="admin-run-requests-0").elements.pop()
    assert switch.value is True
    switch.set_value(False)
    await wait_for(
        lambda: (
            "run_requests"
            not in yaml.safe_load(config_file.read_text())["data_directories"][0]
        )
    )
