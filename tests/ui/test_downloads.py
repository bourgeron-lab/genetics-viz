"""The Download buttons of the cohort and family tabs, in the user simulation."""

import asyncio
import csv
import io
import zipfile

from nicegui import ui
from nicegui.testing import User
from ui_helpers import capture_tables, login, open_tab


async def _download(user: User, mark: str, entry: str = "") -> bytes:
    """Click a Download button (or one of its menu entries) and return the file."""
    target = f"{mark}-{entry}" if entry else mark
    for _ in range(100):
        try:  # find() asserts until the tab has rendered it
            user.find(marker=target)
            break
        except AssertionError:
            await asyncio.sleep(0.05)
    user.find(marker=target).click()
    response = await user.download.next(timeout=5)
    return response.content


def _rows(content: bytes):
    reader = csv.DictReader(io.StringIO(content.decode()), delimiter="\t")
    return reader.fieldnames, list(reader)


async def _open_cohort(user: User, username: str, monkeypatch=None):
    await login(user, username)
    tables = capture_tables(monkeypatch) if monkeypatch else None
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    return tables


async def test_cohort_tabs_download_what_they_show(user):
    await _open_cohort(user, "curator")

    open_tab(user, "cohort-panel-tabs", "Ancestry")
    columns, rows = _rows(await _download(user, "download-ancestry"))
    assert len(rows) == 6
    assert {"region", "PC1", "PC20", "snp_coverage", "admix_EUR"} <= set(columns)

    open_tab(user, "cohort-panel-tabs", "Polygenic Scores")
    columns, rows = _rows(await _download(user, "download-pgs"))
    assert columns[:5] == ["FID", "IID", "sex", "phenotype", "region"]
    assert "height_test" in columns and len(rows) == 6

    open_tab(user, "cohort-panel-tabs", "ROH")
    columns, rows = _rows(await _download(user, "download-roh"))
    assert len(rows) == 6 and "parents_kinship" in columns

    open_tab(user, "cohort-panel-tabs", "STR")
    await user.should_see("flagged genotypes", retries=30)
    columns, rows = _rows(await _download(user, "download-str"))
    assert {(r["IID"], r["locus"]) for r in rows} == {
        ("S3", "HTT"),
        ("S2", "FMR1"),
        ("S6", "ARX_2"),
    }

    open_tab(user, "cohort-panel-tabs", "Sample QC")
    await user.should_see("5 of 6 shown individuals checked", retries=30)
    columns, rows = _rows(await _download(user, "download-sampleqc", "flags"))
    assert "evidence" in columns and "ack_status" in columns
    cross = next(
        r for r in rows if r["code"] == "CRYPTIC_RELATED" and r["sample_a"] == "S2"
    )
    assert cross["ack_note"] == "Known first cousins: their mothers are sisters"
    archive = zipfile.ZipFile(
        io.BytesIO(await _download(user, "download-sampleqc", "zip"))
    )
    names = archive.namelist()
    assert len(names) == 3
    assert {n.split(".")[2] for n in names} == {"flags", "samples", "pairs"}


async def test_cohort_downloads_follow_the_filter_and_the_role(user, monkeypatch):
    tables = await _open_cohort(user, "reader", monkeypatch)
    with user.client:
        tables[0]["on_filter"]({"filters": {"FID": "AU"}})

    open_tab(user, "cohort-panel-tabs", "ROH")
    columns, rows = _rows(await _download(user, "download-roh"))
    assert {r["IID"] for r in rows} == {"S5", "S6"}
    assert "parents_kinship" not in columns and "froh_category" not in columns

    open_tab(user, "cohort-panel-tabs", "Sample QC")
    await user.should_see("2 of 2 shown individuals checked", retries=30)
    columns, rows = _rows(await _download(user, "download-sampleqc", "samples"))
    assert {r["IID"]: r["karyotype"] for r in rows} == {"S5": "Atypical", "S6": "XY"}
    columns, rows = _rows(await _download(user, "download-sampleqc", "pairs"))
    assert "observed_relationship" not in columns
    assert {(r["sample_a"], r["sample_b"]) for r in rows} == {("S5", "S6")}


async def test_family_tabs_download_what_they_show(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)

    columns, rows = _rows(await _download(user, "download-wombat"))
    assert len(rows) == 6 and "actions" not in columns

    open_tab(user, "family-tabs", "SVs")
    columns, rows = _rows(await _download(user, "download-svs"))
    assert len(rows) == 1

    open_tab(user, "family-tabs", "Ancestry")
    columns, rows = _rows(await _download(user, "download-ancestry"))
    assert [r["IID"] for r in rows] == ["S1", "S2", "S3"]
    assert "admix_CSA" in columns

    open_tab(user, "family-tabs", "Polygenic Scores")
    columns, rows = _rows(await _download(user, "download-pgs"))
    # One row per member, as the cohort file, not the table's transpose.
    assert [r["IID"] for r in rows] == ["S1", "S2", "S3"]

    open_tab(user, "family-tabs", "ROH")
    await user.should_see("Compared with the reference", retries=30)
    archive = zipfile.ZipFile(io.BytesIO(await _download(user, "download-roh", "zip")))
    tables = {n.split(".")[2]: n for n in archive.namelist()}
    assert set(tables) == {"summary", "segments", "genes", "kinship"}
    columns, rows = _rows(archive.read(tables["kinship"]))
    # The parents' pair is not shown to readers.
    assert {(r["IID1"], r["IID2"]) for r in rows} == {("S1", "S3"), ("S2", "S3")}

    open_tab(user, "family-tabs", "STR")
    await user.should_see("Cohort DEMO", retries=30)
    columns, rows = _rows(await _download(user, "download-str"))
    assert {r["IID"] for r in rows} == {"S1", "S2", "S3"}
    assert "cohort_percentile" in columns and "validation" in columns

    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    await user.should_not_see("Reading the cohort's tables", retries=30)
    columns, rows = _rows(await _download(user, "download-sampleqc", "flags"))
    assert "evidence" not in columns
    assert {r["code"] for r in rows} == {
        "NOT_PARENT",
        "CONSANGUINITY",
        "POSSIBLE_CONTAMINATION",
        "CRYPTIC_RELATED",
    }
    cross = next(
        r for r in rows if r["code"] == "CRYPTIC_RELATED" and r["sample_a"] == "S2"
    )
    assert cross["ack_status"] == "explained" and cross["ack_note"] == "(hidden)"


async def test_an_empty_export_downloads_nothing(user, monkeypatch):
    tables = await _open_cohort(user, "curator", monkeypatch)
    with user.client:
        tables[0]["on_filter"]({"filters": {"FID": "1002"}})
    open_tab(user, "cohort-panel-tabs", "STR")
    await user.should_see("No flagged genotype among the individuals shown", retries=30)
    user.find(marker="download-str").click()
    await user.should_see("Nothing to download", retries=30)
    assert not user.download.http_responses
    assert user.find(kind=ui.button, marker="download-str").elements
