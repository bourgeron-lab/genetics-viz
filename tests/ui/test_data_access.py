"""Data directories granted per user: the route, the pages and the gene index."""

import asyncio
import sys
from pathlib import Path

import pytest
import yaml
from nicegui import ui
from ui_helpers import capture_tables, login

import datagen

#: DEMO's pedigree in the first directory, EXTRA's in the second.
DEMO_PEDIGREE = "/data-0/cohorts/DEMO/DEMO.pedigree.tsv"
EXTRA_PEDIGREE = "/data-1/cohorts/EXTRA/EXTRA.pedigree.tsv"


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    data = datagen.build_data_dir(tmp_path)
    datagen.add_gene_case(data)
    return data


@pytest.fixture
def config_file(config_file: Path, tmp_path: Path) -> Path:
    """Two directories; the curator may open the second one only."""
    config = yaml.safe_load(config_file.read_text())
    second = datagen.build_second_data_dir(tmp_path)
    config["data_directories"].append({"path": str(second), "description": "second"})
    for entry in config["user_list"]:
        if entry["username"] == "curator":
            entry["data_directories"] = [str(second)]
    config_file.write_text(yaml.safe_dump(config))
    return config_file


async def test_the_route_serves_granted_directories_only(user):
    await login(user, "curator")
    assert (await user.http_client.get(EXTRA_PEDIGREE)).status_code == 200
    assert (await user.http_client.get(DEMO_PEDIGREE)).status_code == 403


async def test_no_grant_means_every_directory(user):
    await login(user, "admin")
    assert (await user.http_client.get(DEMO_PEDIGREE)).status_code == 200
    assert (await user.http_client.get(EXTRA_PEDIGREE)).status_code == 200


async def test_pages_open_on_a_granted_directory(user):
    await login(user, "curator")
    await user.should_see("EXTRA")
    await user.should_not_see("DEMO")
    # A link into the other directory does not switch to it
    await user.open("/goto/0?next=/")
    await user.should_see("EXTRA", retries=30)
    await user.should_not_see("DEMO")


async def test_gene_page_leaves_other_directories_out(user, monkeypatch):
    await login(user, "curator")
    gene_index = sys.modules["genetics_viz.utils.gene_index"]
    await asyncio.to_thread(gene_index.refresh_all, wait=True)
    tables = capture_tables(monkeypatch)
    await user.open("/gene/SHANK3")
    await user.should_see("1 small variant in 1 individual", retries=50)
    cohorts = {r.get("Cohort") for t in tables for r in t["rows"] if r.get("Cohort")}
    assert cohorts == {"EXTRA"}


async def test_no_directory_left_means_no_login(user, config_file):
    config = sys.modules["genetics_viz.config_model"].get_config()
    curator = next(u for u in config.user_list if u.username == "curator")
    curator.data_directories = ["/nowhere"]
    await user.open("/login")
    user.find("Username").type("curator")
    user.find("Password").type("curator-pw")
    user.find(kind=ui.button, content="Sign in").click()
    await user.should_see("No data directory is open to this account", retries=30)
