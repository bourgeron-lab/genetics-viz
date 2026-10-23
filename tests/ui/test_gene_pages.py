"""The gene and variant pages, the header lookup and the data directory switch.

Every test here runs with cohort OTHER beside DEMO and a second data
directory, so a gene's rows span cohorts and directories.
"""

import asyncio
import sys
from pathlib import Path

import pytest
import yaml
from nicegui import ui
from ui_helpers import capture_tables, login, open_tab

import datagen

SHANK3_DE_NOVO = "chr22:50674000:G:A"


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    data = datagen.build_data_dir(tmp_path)
    datagen.add_gene_case(data)
    return data


@pytest.fixture
def config_file(config_file: Path, tmp_path: Path) -> Path:
    config = yaml.safe_load(config_file.read_text())
    second = datagen.build_second_data_dir(tmp_path)
    config["data_directories"].append({"path": str(second), "description": "second"})
    config_file.write_text(yaml.safe_dump(config))
    return config_file


async def _indexed() -> None:
    """Wait for the index the app builds at startup (or build it now)."""
    gene_index = sys.modules["genetics_viz.utils.gene_index"]
    await asyncio.to_thread(gene_index.refresh_all, wait=True)


async def _snapshot_loaded() -> None:
    genetrek = sys.modules["genetics_viz.utils.genetrek"]
    for _ in range(200):
        if genetrek.get_snapshot().ok:
            return
        await asyncio.sleep(0.05)
    raise AssertionError(genetrek.get_snapshot().status_text())


def _table(tables, *columns):
    return next(t for t in tables if all(c in t["columns"] for c in columns))


async def test_gene_page_gathers_every_cohort_and_directory(user, monkeypatch):
    await login(user, "reader")
    await _indexed()
    tables = capture_tables(monkeypatch)
    await user.open("/gene/SHANK3")
    await user.should_see(
        "2 small variants in 3 individuals, 3 families, 3 cohorts, 2 data directories",
        retries=50,
    )
    await user.should_see("cohort tables indexed")

    variants = _table(tables, "Carriers", "Validations")["rows"]
    shared = next(r for r in variants if r["Variant"] == SHANK3_DE_NOVO)
    assert (shared["Carriers"], shared["Cohorts"]) == (3, "DEMO, EXTRA, OTHER")
    assert shared["Validations"] == "absent 1 · present 1 · TODO 1"

    carriers = _table(tables, "Configs", "Validation")["rows"]
    extra = next(r for r in carriers if r["Cohort"] == "EXTRA")
    assert (extra["Data dir"], extra["_dir"]) == ("data2 — second", "1")


async def test_gene_page_tabs_and_link(user, monkeypatch):
    await login(user, "reader")
    await _indexed()
    tables = capture_tables(monkeypatch)
    await user.open("/gene/SCN2A?tab=diag")
    await user.should_see("Diagnostics · 1", retries=50)
    await user.should_see("SVs · 1")
    diagnostics = _table(tables, "Impact", "Users")["rows"]
    assert [(r["Type"], r["Diagnostic"]) for r in diagnostics] == [("SV", "uncertain")]
    open_tab(user, "gene-tabs", "SVs")
    sv = _table(tables, "Locus", "Exonic")["rows"]
    assert [(r["Sample"], r["Validation"]) for r in sv] == [("S3", "")]


async def test_a_previous_symbol_opens_the_current_gene(user, monkeypatch):
    await login(user, "reader")
    await _snapshot_loaded()
    await _indexed()
    tables = capture_tables(monkeypatch)
    await user.open("/gene/KIAA1109")
    await user.should_see("KIAA1109 is now BLTP1; both names are searched.")
    await user.should_see("2 small variants in 1 individual", retries=50)
    variants = _table(tables, "Carriers", "Validations")["rows"]
    assert sorted(r["Genes"] for r in variants) == ["BLTP1", "KIAA1109"]


async def test_str_and_roh_tabs(user, monkeypatch):
    await login(user, "reader")
    await _indexed()
    tables = capture_tables(monkeypatch)
    await user.open("/gene/HTT?tab=str")
    await user.should_see(
        "1 genotype past normal, of 5 individuals genotyped.", retries=50
    )
    rows = _table(tables, "Alleles", "Disease")["rows"]
    assert [(r["Sample"], r["Alleles"]) for r in rows] == [("S3", "18/40")]

    await user.open("/gene/GENEA?tab=roh")
    await user.should_see("ROH · 2", retries=50)


async def test_variant_page_lists_carriers_and_curation_history(user, monkeypatch):
    await login(user, "reader")
    await _indexed()
    tables = capture_tables(monkeypatch)
    await user.open("/variant/22-50674000-G-A")
    await user.should_see(SHANK3_DE_NOVO)
    await user.should_see(
        "Seen in 3 individuals, 3 families, 3 cohorts, 2 data directories.",
        retries=50,
    )
    history = _table(tables, "Ignored", "Comment")["rows"]
    assert [(r["Sample"], r["Validation"], r["Ignored"]) for r in history] == [
        ("S9", "absent", ""),
        ("S9", "present", "yes"),
        ("S3", "present", ""),
    ]
    diagnostics = _table(tables, "Impact", "Users")["rows"]
    assert [r["Diagnostic"] for r in diagnostics] == ["pathogenic"]


async def test_a_variant_seen_nowhere_still_shows_its_curation(user):
    await login(user, "reader")
    await _indexed()
    await user.open("/variant/chr1:1:A:T")
    await user.should_see("Not seen in any cohort Wombat table.", retries=50)
    await user.should_see("No validation recorded for this variant.")
    await user.open("/variant/not-a-variant")
    await user.should_see("Not a variant: not-a-variant")


async def test_header_lookup_opens_the_gene_or_the_variant(user):
    await login(user, "reader")
    user.find(marker="gene-lookup").type("shank3").trigger("keydown.enter")
    await user.should_see("🧬 SHANK3", retries=30)
    user.find(marker="gene-lookup").type("22-50674000-G-A").trigger("keydown.enter")
    await user.should_see(f"🔎 {SHANK3_DE_NOVO}", retries=30)


async def test_goto_switches_the_data_directory(user):
    await login(user, "reader")
    await user.open("/goto/1?next=/cohort/EXTRA")
    await user.should_see("EXTRA", retries=30)
    await user.should_not_see("not found")
    # Only a path on this site is followed
    await user.open("/goto/0?next=//elsewhere.test/")
    await user.should_see("Welcome to Genetics-Viz", retries=30)


async def test_pages_need_a_login(user):
    await user.open("/gene/SHANK3")
    await user.should_see("Sign in", kind=ui.button, retries=30)
