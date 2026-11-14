"""The genome QC views (ghfc-ngs genomeqc) in the Sample QC tabs and headers."""

import asyncio
import csv
import io

import datagen
from nicegui import ui
from nicegui.testing import User
from ui_helpers import capture_tables, login, open_tab

ACK_ROW = (
    "CELL_LINE\t1001\tS3\t\tconfirmed\tKnown LCL, primary DNA requested\tcurator"
    "\t2026-11-14T10:00:00\n"
)


def _genome_table(tables):
    return [t for t in tables if "seq_class" in t["columns"]][-1]


def _issues(tables):
    return [t for t in tables if "severity" in t["columns"]][-1]


def _texts(rows):
    return " ".join(str(v) for r in rows for v in r.values())


async def _download(user: User, mark: str, entry: str) -> bytes:
    target = f"{mark}-{entry}"
    for _ in range(100):
        try:
            user.find(marker=target)
            break
        except AssertionError:
            await asyncio.sleep(0.05)
    user.find(marker=target).click()
    return (await user.download.next(timeout=5)).content


async def _open_family_qc(user, monkeypatch, username, fid="1001"):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open(f"/cohort/DEMO/family/{fid}")
    await user.should_see(f"Family: {fid}", retries=10)
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("DNA and library", retries=30)
    await user.should_not_see("Reading the genome QC tables", retries=30)
    return tables


async def test_curator_family_tab_and_header(user, monkeypatch, data_dir):
    datagen.add_genomeqc_tables(data_dir)
    tables = await _open_family_qc(user, monkeypatch, "curator")
    # The header: CELL_LINE (high) counts with the somalier flags, and the
    # cell line has a badge of its own.
    await user.should_see("QC: 2 high · 1 medium")
    await user.should_see(marker="cell-line-badge")
    await user.should_see("Cell line: S3")
    # The standing warning, with the EBV count for curators.
    await user.should_see(marker="genome-qc-cell-line")
    await user.should_see("152.300 EBV genomes per cell", retries=10)
    # The flag is in the Issues, with the others.
    assert "DNA from a cell line" in _texts(_issues(tables)["rows"])
    table = _genome_table(tables)
    assert "ebv_copies_per_cell" in table["columns"]
    rows = {r["sample"]: r for r in table["rows"]}
    assert set(rows) == {"S1", "S2", "S3"}
    assert rows["S3"]["ebv_copies_per_cell"] == 152.3
    assert rows["S1"]["seq_class_badges"][0]["label"] == "150 bp, PCR"
    assert rows["S1"]["duplicate_share"] == 12.02
    assert rows["S3"]["flags_badges"][0]["label"] == "DNA from a cell line"
    # samtools stats, for curators, through the data route.
    await user.should_see("S3 ↗")
    link = next(
        e
        for e in user.find(ui.link).elements
        if e.props.get("href", "").endswith("S3.gqc_v1.stats.txt")
    )
    response = await user.http_client.get(link.props["href"])
    assert response.status_code == 200 and "samtools stats" in response.text
    await user.should_see(marker="genome-qc-chart")


async def test_reader_gets_no_ebv(user, monkeypatch, data_dir):
    datagen.add_genomeqc_tables(data_dir)
    tables = await _open_family_qc(user, monkeypatch, "reader", fid="AU-X-1")
    table = _genome_table(tables)
    assert "ebv_copies_per_cell" not in table["columns"]
    assert "ebv_copies_per_cell" not in {
        k for r in table["rows"] for k in r
    }  # fmt: skip
    issues = _texts(_issues(tables)["rows"])
    assert "Atypical DNA signal" in issues and "EBV" not in issues
    await user.should_not_see("S5 ↗")
    content = await _download(user, "download-sampleqc", "genome")
    reader = csv.DictReader(io.StringIO(content.decode()), delimiter="\t")
    assert "ebv_copies_per_cell" not in reader.fieldnames
    assert "chr20_depth" in reader.fieldnames
    rows = {r["IID"]: r for r in reader}
    assert rows["S6"]["flags"] == "Low depth"
    assert rows["S5"]["flags"] == "Atypical DNA signal"


async def test_acknowledged_cell_line_keeps_its_badge(user, data_dir):
    datagen.add_genomeqc_tables(data_dir)
    with open(data_dir / "validations" / "sampleqc.tsv", "a") as fh:
        fh.write(ACK_ROW)
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    # Acknowledged: no longer counted, but the caveat stays.
    await user.should_see("QC: 1 high · 1 medium", retries=10)
    await user.should_see("Cell line: S3")


async def test_family_without_cohort_tables(user, monkeypatch, data_dir):
    datagen.add_genomeqc_tables(data_dir, cohort=None)
    tables = await _open_family_qc(user, monkeypatch, "curator")
    await user.should_see("metrics only", retries=10)
    await user.should_not_see(marker="cell-line-badge")
    rows = {r["sample"]: r for r in _genome_table(tables)["rows"]}
    assert rows["S3"]["ebv_copies_per_cell"] == 152.3
    assert rows["S3"]["seq_class_badges"] == [] and rows["S3"]["flags_badges"] == []


async def test_no_genome_qc_no_section(user, monkeypatch):
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    await user.should_not_see("DNA and library")
    await user.should_not_see(marker="cell-line-badge")


async def test_cohort_tab(user, monkeypatch, data_dir):
    datagen.add_genomeqc_tables(data_dir)
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-tabs", "Sample QC")
    await user.should_see("DNA and library", retries=30)
    await user.should_see("5 of 6 shown individuals measured", retries=30)
    await user.should_see(marker="cohort-genome-qc-cell-lines")
    await user.should_see(marker="cohort-genome-qc-chart")
    table = _genome_table(tables)
    assert {r["sample"] for r in table["rows"]} == {"S1", "S2", "S3", "S5", "S6"}
    # The genome QC flags are counted with the others in the flagged table.
    flagged = [t for t in tables if "issues" in t["columns"]][-1]
    assert "S6" in {r["sample"] for r in flagged["rows"]}
    content = await _download(user, "download-sampleqc", "genome")
    reader = csv.DictReader(io.StringIO(content.decode()), delimiter="\t")
    assert "ebv_copies_per_cell" in reader.fieldnames
    assert {r["IID"] for r in reader} == {"S1", "S2", "S3", "S5", "S6"}
