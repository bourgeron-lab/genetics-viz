"""The family report: its sections on the synthetic data, the HTML, the SVGs,
and the curator conclusion file."""

import xml.etree.ElementTree as ET

import datagen
import polars as pl
import pytest

from genetics_viz.models import DataStore
from genetics_viz.utils import report_summaries
from genetics_viz.utils.family_report import (
    build_family_report,
    whole_chromosome_events,
)
from genetics_viz.utils.report_html import render_report_html
from genetics_viz.utils.report_svg import Marker, ordinal, pca_svg, percentile
from genetics_viz.utils.report_svg import pgs_distribution_svg


@pytest.fixture
def store(data_dir, app_config) -> DataStore:
    s = DataStore(data_dir)
    s.load()
    return s


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(report_summaries, "can_write", lambda: True)


def _report(store, fid="1001", cohort=True):
    c = store.cohorts["DEMO"]
    return build_family_report(
        store.data_dir, fid, c.families[fid], c if cohort else None, "curator"
    )


def test_sections_follow_the_requested_order(store):
    report = _report(store)
    assert [s.key for s in report.sections] == [
        "pedigree", "ancestry", "sample_qc", "roh", "pgs", "aneuploidy", "upd",
        "diagnostics",
    ]  # fmt: skip
    assert report.document_title.startswith("DEMO_1001_report_")
    assert report.pipeline.startswith("ghfc-ngs ")


def test_pedigree_and_ancestry(store):
    report = _report(store)
    pedigree = report.section("pedigree").rows
    assert [(r["iid"], r["father"], r["mother"], r["phenotype"]) for r in pedigree] == [
        ("S1", "—", "—", "Unaffected"),
        ("S2", "—", "—", "Unaffected"),
        ("S3", "S1", "S2", "Affected"),
    ]
    assert [r["karyotype"] for r in pedigree] == ["XY", "XX", "XY"]
    ancestry = report.section("ancestry")
    assert {r["region"] for r in ancestry.rows} == {"EUR"}
    assert ancestry.rows[0]["admixture"] == "EUR 80% · CSA 20%"
    assert ancestry.glance == "EUR (3 individuals)"
    assert len(ancestry.plots) == 2  # the whole panel, and zoomed


def test_sample_qc_uses_curator_wording_and_keeps_reviews(store):
    section = _report(store).section("sample_qc")
    issues = [r["issue"] for r in section.rows if not r.get("_muted")]
    assert "Not the declared parent" in issues
    assert "Possible contamination" in issues  # the withdrawn ack does not count
    reviewed = [r for r in section.rows if r.get("_muted")]
    assert reviewed and all(
        r["review"].startswith(("Explained", "Confirmed")) for r in reviewed
    )
    assert section.attention


def test_without_cohort_cross_family_flags_are_not_checked(store):
    section = _report(store, cohort=False).section("sample_qc")
    assert "Relatedness with other families was not checked" in section.note
    assert all(" – S5" not in r["samples"] for r in section.rows)


def test_roh_reads_the_grade_and_relationship(store):
    rows = {r["iid"]: r for r in _report(store).section("roh").rows}
    assert rows["S3"]["grade"] == "High"
    assert "second-degree-like" in rows["S3"]["_grade_note"]
    assert rows["S3"]["pct"].startswith("96th (vs EUR")
    assert rows["S1"]["_grade_note"] == ""


def test_pgs_missing_then_present(store, data_dir):
    section = _report(store).section("pgs")
    assert not section.available
    datagen.add_report_case(data_dir)
    section = _report(store).section("pgs")
    assert section.available and len(section.plots) == 4
    rows = {r["iid"]: r for r in section.rows}
    assert rows["S3"]["autism_grove2026"] == "z +2.40 · >99th"
    assert rows["S3"]["adhd_demontis2023"] == "—"
    assert rows["S1"]["int_savage2018"] == "z +1.10 · 86th"
    assert section.glance == "S3 Autism >99th"


def test_aneuploidy_sex_chromosomes_and_whole_chromosome_gain(store, data_dir):
    au = _report(store, "AU-X-1").section("aneuploidy")
    assert au.glance == "S5: XXY; S6: XY, declared female"
    assert _report(store).section("aneuploidy").note.startswith("No whole-chromosome")
    datagen.add_report_case(data_dir)
    section = _report(store).section("aneuploidy")
    rows = {r["iid"]: r for r in section.rows}
    assert rows["S3"]["autosomes"] == "+21 (100%)"
    assert rows["S1"]["autosomes"] == "none detected"
    assert section.attention


def test_whole_chromosome_events_need_most_of_the_chromosome():
    calls = pl.DataFrame(
        {
            "chr": ["chr18", "chr18", "chr2"],
            "start": [0.0, 40e6, 1e6],
            "end": [30e6, 80e6, 2e6],
            "ratio": [-0.9, -0.8, 0.5],
            "sample": ["A", "A", "A"],
        }
    )
    assert whole_chromosome_events(calls, 0.75) == {"A": ["−18 (87%)"]}  # 70 of 80.4 Mb
    assert whole_chromosome_events(calls, 0.95) == {}


def test_upd_candidate_from_both_signals(store, data_dir):
    assert _report(store).section("upd").rows == []
    datagen.add_upd_case(data_dir)
    section = _report(store).section("upd")
    assert [
        (r["iid"], r["chrom"], r["evidence"], r["validation"]) for r in section.rows
    ] == [("S3", "chr2", "ROH + Sample QC", "not reviewed")]
    assert section.attention


def test_diagnostics_keep_pathogenic_and_uncertain(store, data_dir):
    datagen.add_report_case(data_dir)
    section = _report(store).section("diagnostics")
    assert [(r["verdict"], r["kind"]) for r in section.rows] == [
        ("Pathogenic", "SNV"),
        ("Uncertain significance", "SV"),
    ]
    snv, sv = section.rows
    assert snv["validation"] == "S3: present (de novo)"
    assert sv["finding"].startswith("DEL chr2:165114001-165316000 (2q")
    assert sv["genes"] == "SCN2A"
    assert section.glance == "1 pathogenic · 1 uncertain significance"


def test_family_without_results_degrades(store):
    report = _report(store, "1002")
    by_key = {s.key: s for s in report.sections}
    assert not by_key["sample_qc"].available
    assert not by_key["aneuploidy"].available
    assert by_key["roh"].available
    titles = [t for t, _ in report.methods]
    assert "Sample quality control" not in titles and "Genetic ancestry" in titles


def test_methods_fill_the_run_values(store):
    text = dict(_report(store).methods)["Genetic ancestry"]
    assert "bundle 1.1.0, genotype filters dp10gq20" in text
    cnv = dict(_report(store).methods)["Copy-number variants and aneuploidy"]
    assert "n/a" not in cnv and "fixed-size bins" in cnv


def test_html_is_escaped_and_svgs_parse(store, data_dir, curator):
    datagen.add_report_case(data_dir)
    report_summaries.save_family_summary(
        data_dir, "1001", "Line one <script>x</script>\nLine two", user="curator"
    )
    report = _report(store)
    html = render_report_html(report)
    assert html.startswith("<!DOCTYPE html>")
    assert "<script>x</script>" not in html
    assert "&lt;script&gt;x&lt;/script&gt;\nLine two" in html
    assert f"<title>{report.document_title}</title>" in html
    for section in report.sections:
        for plot in section.plots:
            ET.fromstring(plot.svg)
    assert html.count("<svg") == 6


def test_svg_builders_escape_labels():
    svg = pgs_distribution_svg([Marker("<b>", 5.0, "#000"), Marker("S2", -0.5, "#111")])
    root = ET.fromstring(svg)
    assert any(
        "<b>" in (t.text or "") for t in root.iter("{http://www.w3.org/2000/svg}text")
    )
    ET.fromstring(pca_svg([("S1", 0.1, 0.2, "#000")], [(0.0, 0.0, "EUR")], None))
    assert pca_svg([], None, None) == ""


def test_percentile_wording():
    assert ordinal(percentile(0.5)) == "69th"
    assert ordinal(percentile(-3)) == "<1st"
    assert ordinal(percentile(3)) == ">99th"
    assert ordinal(percentile(1.0)) == "84th"


def test_summaries_role_guard_and_last_row_wins(tmp_path, curator, monkeypatch):
    report_summaries.save_family_summary(tmp_path, "F1", "first", user="a")
    report_summaries.save_family_summary(tmp_path, "F1", "a\\b\nc", user="b")
    latest = report_summaries.latest_family_summary(tmp_path, "F1")
    assert (latest.text, latest.user) == ("a\\b\nc", "b")
    report_summaries.save_family_summary(tmp_path, "F1", "  ", user="b")
    assert report_summaries.latest_family_summary(tmp_path, "F1") is None
    monkeypatch.setattr(report_summaries, "can_write", lambda: False)
    with pytest.raises(PermissionError):
        report_summaries.save_family_summary(tmp_path, "F1", "x", user="c")
