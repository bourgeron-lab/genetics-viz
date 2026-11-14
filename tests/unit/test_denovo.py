"""The denovo step's tables (utils/denovo.py) and their review count."""

import os

import datagen
import polars as pl
import pytest

from genetics_viz.utils import denovo
from genetics_viz.utils import review_counts as rc


@pytest.fixture
def demo(data_dir):
    """1001's denovo tables, and the DEMO cohort (name and path are enough)."""
    datagen.add_denovo_calls(data_dir)

    class Demo:
        name = "DEMO"
        path = data_dir / "cohorts" / "DEMO"

    return Demo()


def _explain_not_parent(data_dir):
    with open(data_dir / "validations" / "sampleqc.tsv", "a") as f:
        f.write(
            "\t".join(["NOT_PARENT", "1001", "S1", "S3", "explained", "Adopted",
                       "curator", "2026-11-10T10:00:00"]) + "\n"
        )  # fmt: skip


def _rows(data):
    """1001's Wombat de novo rows, as the tab builds them, with evidence."""
    return [
        {"Variant": f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}", "sample": "S3"}
        for r in data.rows()
    ]


def test_family_reads_the_cohort_table_when_as_new(data_dir, demo):
    data = denovo.load_family_denovo(data_dir, "1001", demo)
    assert data.level == "cohort" and data.has_recurrence
    assert data.label == "dnm_v2" and data.wombat_config == "rare_denovo"
    assert data.qc_of("S3")["qc"] == "low"
    # Without the cohort, or with an older cohort table: the family's own
    alone = denovo.load_family_denovo(data_dir, "1001")
    assert alone.level == "family" and not alone.has_recurrence
    old = os.path.getmtime(datagen.family_dir(data_dir, "1001") / "denovo") - 100
    for p in (demo.path / "denovo").iterdir():
        os.utime(p, (old, old))
    assert denovo.load_family_denovo(data_dir, "1001", demo).level == "family"


def test_family_id_is_read_as_a_string(data_dir, demo):
    data = denovo.load_cohort_denovo(demo)
    assert data.dnm["family_id"].dtype == pl.Utf8
    assert set(data.dnm["family_id"].to_list()) == {"1001"}


def test_evidence_and_strict_filter(data_dir, demo):
    data = denovo.load_family_denovo(data_dir, "1001", demo)
    rows = _rows(data) + [{"Variant": "chr9:1:A:C", "sample": "S3"}]
    assert denovo.add_evidence(rows, data) == 1
    shank3, plain, unknown = rows
    assert shank3["dnm_strict"] == "yes" and shank3["dnm_origin"] == "paternal"
    assert shank3["dnm_child_reads"] == "14/29" and shank3["_dnm"]["genes"] == "SHANK3"
    assert plain["dnm_flags"] == "PARENT_ALT_POSSIBLE,TR_HOMOPOLYMER"
    assert [b["label"] for b in plain["DnmFlagBadges"]] == [
        "Parent alt 1",
        "Tandem repeat",
    ]
    assert "Wright 2019" in plain["DnmFlagBadges"][0]["tooltip"]
    assert unknown["dnm_flags"] == denovo.NOT_ASSESSED

    f = denovo.DenovoFilter()
    # Strict: the flagged one goes, the one not assessed stays
    assert denovo.apply_denovo_filter(rows, f) == [shank3, unknown]
    f.mode["value"] = "all"
    assert len(denovo.apply_denovo_filter(rows, f)) == 3
    f.hide["value"] = ["TR_HOMOPOLYMER"]
    assert plain not in denovo.apply_denovo_filter(rows, f)
    f.hide["value"], f.only["value"] = [], ["PARENT_ALT_POSSIBLE"]
    assert denovo.apply_denovo_filter(rows, f) == [plain]


def test_hidden_curated_calls(data_dir, demo):
    data = denovo.load_family_denovo(data_dir, "1001", demo)
    rows = _rows(data)
    denovo.add_evidence(rows, data)
    rows[1].update(Validation="present", ValidationInheritance="de novo")
    shown = denovo.apply_denovo_filter(rows, denovo.DenovoFilter())
    assert denovo.hidden_curated(rows, shown) == [rows[1]]
    # A call validated as inherited is not a lost de novo call
    rows[1]["ValidationInheritance"] = "maternal"
    assert denovo.hidden_curated(rows, shown) == []


def test_curation_check_finds_hidden_and_missing(data_dir, demo):
    with open(data_dir / "validations" / "snvs.tsv", "a") as f:
        f.write("\t".join(["1001", "chr5:5000:A:G", "S3", "curator", "de novo",
                           "present", "", "0", "2026-11-14T10:00:00"]) + "\n")  # fmt: skip
    data = denovo.load_cohort_denovo(demo)
    curation = denovo.load_curation(data_dir, {"1001"})
    issues = {
        (r["Variant"], r["problem"])
        for r in denovo.curation_check(data.rows(), curation, {"1001"})
    }
    assert issues == {("chr2:2000:G:T", "hidden"), ("chr5:5000:A:G", "missing")}
    table = denovo.strict_vs_curation(data.rows(), curation)
    assert table["no"]["present"] == 1
    flags = {r["flag"]: r for r in denovo.flag_summary(data.rows(), curation)}
    assert flags["TR_HOMOPOLYMER"]["present"] == 1 and flags["TR_HOMOPOLYMER"]["n"] == 1


def test_child_rows_and_summary(data_dir, demo):
    data = denovo.load_cohort_denovo(demo)
    (child,) = denovo.child_rows(data)
    assert (child["FID"], child["child"], child["qc"]) == ("1001", "S3", "low")
    assert (child["n_candidates"], child["n_strict"], child["n_strict_snv"]) == (
        2,
        1,
        1,
    )
    assert child["paternal_share"] == 1.0 and child["ti_fraction"] == 1.0
    assert child["snv_per_gb"] == round(1 / 2.3, 2)
    assert denovo.child_rows(data, {("1001", "S1")}) == []
    summary = denovo.cohort_summary([child])
    assert summary["n_strict"] == 1 and summary["implied_snv_mu"] > 0
    assert {r["metric"] for r in denovo.REFERENCES} <= set(summary) | {"ti_fraction"}


def test_source_is_the_steps_config(data_dir, demo):
    from genetics_viz.utils.wombat import find_wombat_files

    files = find_wombat_files(datagen.family_dir(data_dir, "1001") / "wombat", "1001")
    data = denovo.load_family_denovo(data_dir, "1001", demo)
    assert denovo.denovo_source(files, data)["wombat_config"] == "rare_denovo"
    assert denovo.denovo_source(files, None)["wombat_config"] == "rare_denovo"
    assert (
        denovo.denovo_source([f for f in files if "denovo" not in str(f)], None) is None
    )


def test_review_counts_strict_candidates_only(data_dir, demo):
    _explain_not_parent(data_dir)
    count = rc.denovo_to_review(data_dir, "1001", demo)
    # SHANK3 (strict, missense) counts; PLAIN is flagged (and validated)
    assert count.n == 1 and "strict de novo candidate" in count.tooltip
    assert "LCR" not in count.tooltip


def test_cohort_status_and_pills(data_dir, demo):
    verdicts, issues = rc.denovo_cohort_status(demo, data_dir, with_curation=True)
    assert verdicts == {("1001", "S3"): "low"}
    assert issues == [("1001", "S3")]
    pills = rc.denovo_cohort_pills(verdicts, issues, {("1001", "S3")})
    assert [p.text for p in pills] == ["1", "1"]
    assert rc.denovo_cohort_pills(verdicts, issues, {("1002", "S4")}) == []
    _, none = rc.denovo_cohort_status(demo, data_dir, with_curation=False)
    assert none == []


def test_child_qc_alerts(data_dir, demo):
    rows = denovo.load_child_qc(data_dir, "1001", demo)
    assert denovo.qc_alerts(rows) == [("S3", "low")]
    # The family table has no verdict: no alert
    assert denovo.qc_alerts(denovo.load_child_qc(data_dir, "1001")) == []


# ---------------------------------------------------------------------------
# With the genome QC (ghfc-ngs genomeqc): classes, cell lines, germline checks
# ---------------------------------------------------------------------------


@pytest.fixture
def demo_gqc(data_dir):
    """As ``demo``, the cohort QC table written with a genome QC: S3 a cell line."""
    datagen.add_denovo_calls(data_dir, genome_qc=True)

    class Demo:
        name = "DEMO"
        path = data_dir / "cohorts" / "DEMO"

    return Demo()


def test_child_rows_carry_the_germline_checks(data_dir, demo_gqc):
    data = denovo.load_cohort_denovo(demo_gqc)
    (child,) = denovo.child_rows(data)
    assert child["qc"] == "high" and child["cell_line"]
    assert child["seq_class"] == "rl150"
    assert child["SeqClassBadges"][0]["label"] == "150 bp"
    assert [b["label"] for b in child["QcFlagBadges"]] == [
        "Cell line",
        "Not germline-like",
    ]
    assert child["paternal_p"] == pytest.approx(1.2e-42)
    assert child["expected_germline_share"] == 0.0424
    assert child["n_phased_strict_snv"] == 522
    # The app's own paternal share (all strict, recounted) stays apart from
    # the pipeline's (strict SNVs only).
    assert child["paternal_share"] == 1.0 and child["paternal_share_snv"] == 0.4713


def test_germline_wording(data_dir, demo_gqc):
    (q,) = denovo.load_child_qc(data_dir, "1001", demo_gqc)
    assert denovo.is_cell_line(q)
    assert denovo.germline_text(q) == (
        "expected germline share 4%; paternal 47% of 522 phased strict SNVs"
        " (p = 1.2e-42); spectrum cosine 0.88"
    )
    assert denovo.qc_reference_text("cohort:rl150") == (
        "the cohort's 150 bp children (median ± 4 scaled MADs, cell lines left out)"
    )
    assert denovo.qc_reference_text("fixed").startswith("the fixed range")
    assert denovo.qc_flag_badges("NEW_FLAG")[0]["label"] == "NEW_FLAG"
    assert denovo.germline_text({}) == ""
    # Within the range, no expected germline share (it would read as somatic)
    ok = {**q, "qc": "ok", "qc_flags": None}
    assert denovo.germline_text(ok).startswith("paternal 47%")


def test_a_flag_alone_raises_the_alert():
    # A polyclonal cell line can have a normal count: still an alert.
    q = {"child": "S3", "qc": "ok", "qc_flags": "CELL_LINE"}
    assert denovo.qc_alerts([q]) == [("S3", "ok")]
    assert denovo.qc_alerts([{**q, "qc_flags": None}]) == []
    assert denovo.is_cell_line({"genomeqc_flags": "EBV_TRACE,CELL_LINE"})
