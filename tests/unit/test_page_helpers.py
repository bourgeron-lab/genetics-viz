"""Data helpers behind the cohort and family tabs, and the change monitor."""

import asyncio
import os

import polars as pl

from genetics_viz.config_model import DataDirectoryConfig
from genetics_viz.pages.cohort.components import (
    ancestry_tab,
    cohort_ancestry_tab,
    cohort_pgs_tab,
    roh_tab,
    stats_panel,
)
from genetics_viz.utils import change_monitor, data
from genetics_viz.utils.ancestry import find_ancestry_files, find_cohort_ancestry_files
from genetics_viz.utils.roh import find_roh_files


class DemoCohort:
    name = "DEMO"

    def __init__(self, data_dir):
        self.path = data_dir / "cohorts" / "DEMO"


def test_family_ancestry_rows_join_predictions(data_dir):
    rows = ancestry_tab._load_frames(find_ancestry_files(data_dir, "1001"))
    assert [r["IID"] for r in rows] == ["S1", "S2", "S3"]
    assert rows[0]["region"] == "EUR" and rows[0]["PC1"] == 0.01


def test_cohort_ancestry_and_pgs_survive_text_family_ids(data_dir):
    files = find_cohort_ancestry_files(DemoCohort(data_dir))
    rows = cohort_ancestry_tab._load_rows(files)
    assert {r["family_id"] for r in rows} == {"1001", "1002", "AU-X-1"}
    assert dict(cohort_ancestry_tab._region_counts(rows)) == {"EUR": 6}
    by_sample, traits, regions = cohort_pgs_tab._load(files)
    assert set(by_sample) == {"S1", "S2", "S3", "S4", "S5", "S6"}
    assert traits == ["bmi_test", "height_test"] and regions == ["EUR"]
    assert by_sample["S1"]["bmi_test"] is None  # the literal 'nan' is nulled


def test_roh_tab_loads_optional_tables(data_dir):
    tables = roh_tab._load(find_roh_files(data_dir, "1002"))
    assert tables["froh"].height == 1 and tables["kinship"] is None


def test_stats_counts_unique_variants_per_sample(tmp_path):
    path = tmp_path / "w.tsv"
    pl.DataFrame(
        {
            "#CHROM": ["1", "1", "2"],
            "POS": [1, 1, 2],
            "REF": ["A", "A", "C"],
            "ALT": ["T", "T", "G"],
            "sample": ["S1", "S1", "S2"],
            "VEP_SYMBOL": ["GENEA", "GENEA", "GENEB,GENEC"],
        }
    ).write_csv(path, separator="\t")
    assert stats_panel._count_variants_per_sample(path, None) == {"S1": 1, "S2": 1}
    assert stats_panel._count_variants_per_sample(path, {"GENEC"}) == {"S2": 1}


def test_change_monitor_reports_and_survives_a_failing_subscriber(data_dir):
    data.init_all_data_stores([DataDirectoryConfig(str(data_dir))])
    received = []

    async def good(report):
        received.append(report)

    async def bad(report):
        raise RuntimeError("subscriber failure")

    unsubscribe = change_monitor.subscribe(good)
    change_monitor.subscribe(bad)
    assert asyncio.run(change_monitor.check_now()) == []

    pedigree = data_dir / "cohorts" / "DEMO" / "DEMO.pedigree.tsv"
    st = pedigree.stat()
    os.utime(pedigree, (st.st_atime, st.st_mtime + 10))
    reports = asyncio.run(change_monitor.check_now())
    assert [r.modified_cohorts for r in reports] == [["DEMO"]] and len(received) == 1

    unsubscribe()
    os.utime(pedigree, (st.st_atime, st.st_mtime + 20))
    asyncio.run(change_monitor.check_now())
    assert len(received) == 1


def test_ancestry_rows_keep_the_pcs_order_whatever_the_predictions_order(data_dir):
    """Polars 2 guarantees no join order unless asked; the tabs ask for it."""
    files = find_ancestry_files(data_dir, "1001")
    predictions = files.ancestry_path.read_text().splitlines()
    files.ancestry_path.write_text(
        "\n".join([predictions[0], *reversed(predictions[1:])]) + "\n"
    )
    assert [r["IID"] for r in ancestry_tab._load_frames(files)] == ["S1", "S2", "S3"]
