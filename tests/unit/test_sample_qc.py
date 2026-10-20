"""Sample QC: reading, flag wording, shaping and plots (utils/sample_qc*.py)."""

import json

import datagen
import polars as pl

from genetics_viz.utils import sample_qc as qc
from genetics_viz.utils import sample_qc_plots as plots


class Cohort:
    def __init__(self, data_dir):
        self.name = "DEMO"
        self.path = data_dir / "cohorts" / "DEMO"


def test_discovery(data_dir):
    assert qc.probe_qc_data(data_dir, "1001")
    assert qc.probe_qc_data(data_dir, "AU-X-1")
    assert not qc.probe_qc_data(data_dir, "1002")
    cohort = Cohort(data_dir)
    assert qc.probe_cohort_qc_data(cohort)
    assert qc.find_qc_labels(qc.get_cohort_qc_dir(cohort), "DEMO") == [datagen.QC_LABEL]
    assert (
        qc.cohort_report_url(cohort, datagen.QC_LABEL, "/data-0")
        == "/data-0/cohorts/DEMO/qc/DEMO.somalier_hg38.html"
    )


def test_dotted_label_is_kept_whole(tmp_path):
    (tmp_path / "F1.somalier.hg38.v2.samples.tsv").write_text("family_id\tsample\n")
    assert qc.find_qc_labels(tmp_path, "F1") == ["somalier.hg38.v2"]


def test_read_keeps_ids_as_text_and_adds_missing_columns(tmp_path):
    rows = [f"{88171200000 + i}\tS{i}\t." for i in range(50)] + ["AU-X\tSx\tFLAG"]
    path = tmp_path / "x.samples.tsv"
    # A stub run's cut-down header.
    path.write_text("\n".join(["family_id\tsample\tflags", *rows]) + "\n")
    df = qc.read_qc_table(path, qc.SAMPLE_COLUMNS)
    assert df.schema["family_id"] == pl.String
    assert df["flags"][0] is None and df["flags"][-1] == "FLAG"
    assert set(qc.SAMPLE_COLUMNS) <= set(df.columns)
    assert df["x_copies"].to_list()[0] is None


def test_read_none_and_empty(tmp_path):
    assert qc.read_qc_table(None) is None
    (tmp_path / "e.tsv").write_text("")
    assert qc.read_qc_table(tmp_path / "e.tsv") is None


def test_family_tables(data_dir):
    data = qc.load_family_qc(data_dir, "1001")
    assert data.label == datagen.QC_LABEL
    assert [r["sample"] for r in data.sample_rows()] == ["S1", "S2", "S3"]
    assert [f.code for f in data.flags] == ["NOT_PARENT", "POSSIBLE_CONTAMINATION"]
    not_parent = data.flags[0]
    assert not_parent.samples == ["S1", "S3"] and not not_parent.is_cross_family
    pairs = qc.pair_rows(data.pairs)
    mismatch = next(p for p in pairs if p.status == "mismatch")
    assert (mismatch.sample_a, mismatch.sample_b) == ("S1", "S3")
    assert mismatch.declared and mismatch.po_rate == 0.03


def test_flag_wording():
    karyo = qc.flag_info("KARYOTYPE_XXYY")
    assert karyo.category == "karyotype" and karyo.label == "Karyotype XXYY"
    assert karyo.reader_label == "Sex chromosome anomaly"
    unknown = qc.flag_info("SOMETHING_NEW")
    assert unknown.label == "SOMETHING_NEW" and unknown.category == "sample_quality"
    assert qc.severity_info("bogus")["key"] == "info"
    assert (
        qc.flag_info("NOT_PARENT").reader_label == "Pedigree relationship inconsistent"
    )


def test_reader_wording_hides_evidence():
    f = qc.FlagRow("family", "1001", "S1", "S3", "NOT_PARENT", "high", "S1 is not...")
    reader = qc.flag_badge(f, show_detail=False)
    assert reader["label"] == "Pedigree relationship inconsistent"
    assert "S1 is not" not in reader["tooltip"]
    curator = qc.flag_badge(f, show_detail=True)
    assert curator["label"] == "Not the declared parent"
    assert "Evidence: S1 is not" in curator["tooltip"]


def test_cross_family_flag_keys_split_the_families():
    f = qc.FlagRow("cohort", "1001,AU-X-1", "S2", "S5", "CRYPTIC_RELATED", "high", "")
    assert f.keys == [("1001", "S2"), ("AU-X-1", "S5")]
    assert f.is_cross_family and f.other_family("1001") == "AU-X-1"
    by_key = qc.flags_by_key([f])
    assert set(by_key) == {("1001", "S2"), ("AU-X-1", "S5")}


def test_family_flags_gain_the_cohort_ones_once(data_dir):
    family = qc.load_family_qc(data_dir, "1001")
    cohort = qc.load_cohort_qc(Cohort(data_dir))
    flags = qc.family_flags_with_cohort(family, cohort, "1001")
    assert sorted(f.code for f in flags) == [
        "CRYPTIC_RELATED",
        "NOT_PARENT",
        "POSSIBLE_CONTAMINATION",
    ]
    assert qc.family_flags_with_cohort(family, None, "1001") == family.flags


def test_counts_restrict_to_and_by_individual(data_dir):
    cohort = qc.load_cohort_qc(Cohort(data_dir))
    assert qc.count_by_severity(cohort.flags) == {"high": 4, "medium": 1, "info": 0}
    only_au = {("AU-X-1", "S5"), ("AU-X-1", "S6")}
    assert qc.count_by_severity(cohort.flags, only_au)["high"] == 3
    # S1, S2, S3 (high via NOT_PARENT), S5, S6: S3's worst is high.
    assert qc.worst_by_individual(cohort.flags) == {"high": 5, "medium": 0, "info": 0}
    assert qc.severity_summary({"high": 2, "medium": 1}) == "2 high · 1 medium"
    assert qc.severity_summary({"high": 0, "medium": 0, "info": 3}) == ""


def test_category_badges_are_neutral_and_ordered(data_dir):
    cohort = qc.load_cohort_qc(Cohort(data_dir))
    s5 = qc.flags_by_key(cohort.flags)[("AU-X-1", "S5")]
    badges = qc.category_badges(s5, show_detail=False)
    assert [b["label"] for b in badges] == ["Across families", "Sex chromosomes"]
    assert all("\n" in b["tooltip"] for b in badges)
    assert "observed second-degree" not in badges[0]["tooltip"]


def test_cross_family_graph_counts_relatives():
    pairs = [
        qc.PairRow("A", "a1", "B", "b1", 0.05, 0.15, 1, 0, "unrelated", "third-degree",
                   "mismatch", None),
        qc.PairRow("A", "a1", "C", "c1", 0.05, 0.15, 1, 0, "unrelated", "third-degree",
                   "mismatch", None),
        qc.PairRow("A", "a1", "D", "d1", 0.05, 0.15, 1, 0, "unrelated", "third-degree",
                   "mismatch", None),
        qc.PairRow("A", "a1", "A", "a2", 0.25, 0.0, 1, 0.25, "parent-offspring",
                   "parent-offspring", "ok", None),
    ]  # fmt: skip
    graph = qc.cross_family_graph(pairs)
    assert len(graph.edges) == 3
    assert [(h.sample, h.degree) for h in graph.hubs(3)] == [("a1", 3)]
    restricted = qc.cross_family_graph(pairs, restrict_to={("B", "b1")})
    assert len(restricted.edges) == 1


def test_denovo_risk(data_dir):
    family = qc.load_family_qc(data_dir, "1001")
    risky = qc.denovo_risk(family.flags, ["S1", "S2", "S3"])
    assert {f.code for f in risky} == {"NOT_PARENT", "POSSIBLE_CONTAMINATION"}
    au = qc.load_family_qc(data_dir, "AU-X-1")
    assert qc.denovo_risk(au.flags, ["S5", "S6"]) == []


def test_mendelian_fence():
    rows = [
        {"trio_mendelian_error_rate": v}
        for v in (0.0008, 0.0008, 0.0009, 0.0008, 0.0027)
    ]
    fence = qc.mendelian_fence(rows)
    assert 0.0009 < fence < 0.0027
    assert qc.mendelian_fence(rows[:2]) is None


def test_relatedness_scatter(data_dir):
    pairs = qc.pair_rows(qc.load_family_qc(data_dir, "1001").pairs)
    curator = plots.relatedness_scatter_option(pairs, title="t", show_detail=True)
    names = [s["name"] for s in curator["series"]]
    assert "Parent–offspring" in names and "Mismatch" in names
    mismatch = next(s for s in curator["series"] if s["name"] == "Mismatch")
    assert [d["name"] for d in mismatch["data"]] == ["S1 – S3"]
    lines = [
        m["name"]
        for s in curator["series"]
        for m in s.get("markLine", {}).get("data", [])
    ]
    assert "PO | FS" in lines and any("1st degree" in n for n in lines)
    reader = plots.relatedness_scatter_option(pairs, title="t", show_detail=False)
    assert {s["name"] for s in reader["series"]} == {"Pairs", "Mismatch"}
    assert "observed" not in json.dumps(reader)


def test_kinship_matrix(data_dir):
    pairs = qc.pair_rows(qc.load_family_qc(data_dir, "1001").pairs)
    option = plots.kinship_matrix_option(["S1", "S2", "S3"], pairs, show_detail=True)
    data = option["series"][0]["data"]
    assert len(data) == 6  # three pairs, both triangles
    outlined = [
        d["value"][:2] for d in data if d["itemStyle"]["borderColor"] == "#c62828"
    ]
    assert sorted(outlined) == [[0, 2], [2, 0]]


def test_sex_and_quality_plots(data_dir):
    rows = qc.load_cohort_qc(Cohort(data_dir)).sample_rows()
    points = [
        dict(qc.sample_point(r), flagged=r["sample"] in ("S5", "S6")) for r in rows
    ]
    option = plots.sex_scatter_option(points, show_detail=True)
    areas = option["series"][0]["markArea"]["data"]
    assert {a[0]["name"] for a in areas} >= {"XX", "XY", "XXY", "X0"}
    flagged = next(s for s in option["series"] if s["name"] == "Flagged")
    assert {d["name"] for d in flagged["data"]} == {"S5", "S6"}
    het = plots.x_het_option(points, show_detail=False)
    assert "karyotype: XXY" not in json.dumps(het)
    strip = plots.metric_strip_option(
        [{"name": "S3", "value": 0.071, "flagged": True}],
        background=[0.02, 0.018],
        title="",
        axis_name="CHARR",
        lines=[(0.05, "flag", "#c62828")],
    )
    assert strip["xAxis"]["max"] >= 0.071
    json.dumps(option), json.dumps(strip)


def test_cross_family_graph_option_marks_hubs():
    pairs = [
        qc.PairRow("A", "a1", f"F{i}", f"x{i}", 0.05, 0.15, 1, 0, "unrelated",
                   "third-degree", "mismatch", None)
        for i in range(5)
    ]  # fmt: skip
    option = plots.cross_family_graph_option(
        qc.cross_family_graph(pairs), show_detail=False
    )
    nodes = option["series"][0]["data"]
    hub = next(n for n in nodes if n["name"] == "a1")
    assert hub["label"]["show"] and hub["value"] == 5
    assert "third" not in json.dumps(option["series"][0]["links"])
