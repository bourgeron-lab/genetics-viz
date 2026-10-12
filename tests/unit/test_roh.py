"""Runs of homozygosity: reading, readings and kinship (utils/roh.py)."""

import polars as pl
import pytest

from genetics_viz.utils import roh


def test_read_roh_table_keeps_ids_as_text_and_strips_chr(tmp_path):
    rows = [f"S{i}\tchr{1 + i % 22}\t{1000 + i}" for i in range(150)] + [
        "Sx\tchrM\tAU-ANG-FER-MPX3"
    ]
    path = tmp_path / "x.roh.tsv"
    path.write_text("\n".join(["IID\tchrom\tfamily_id", *rows]) + "\n")
    df = roh.read_roh_table(path)
    assert df.schema["family_id"] == pl.String
    assert df["chrom"][0] == "1" and df["chrom"][-1] == "MT"


def test_read_roh_table_none_and_empty(tmp_path):
    assert roh.read_roh_table(None) is None
    (tmp_path / "e.tsv").write_text("")
    assert roh.read_roh_table(tmp_path / "e.tsv") is None


def test_family_and_cohort_discovery(data_dir):
    files = roh.find_roh_files(data_dir, "1001")
    assert {"froh", "roh", "roh_genes", "kinship"} <= set(files.paths)
    assert roh.find_roh_files(data_dir, "1002").path("kinship") is None  # singleton

    class Cohort:
        name = "DEMO"
        path = data_dir / "cohorts" / "DEMO"

    assert roh.probe_cohort_roh_data(Cohort)
    assert roh.find_cohort_roh_files(Cohort).path("roh_windows") is not None
    assert roh.probe_roh_data(data_dir, "1001") and not roh.probe_roh_data(
        data_dir, "zzz"
    )


def test_froh_qc_pass_is_boolean(data_dir):
    froh = roh.read_roh_table(roh.find_roh_files(data_dir, "1001").path("froh"))
    assert froh.schema["qc_pass"] == pl.Boolean
    assert froh.schema["upd_candidate_chrom"] == pl.String


def test_grades():
    assert [g.threshold for g in roh.threshold_grades()] == [0.03, 0.10, 0.17]
    assert roh.grade_for("first-cousin-like").grade == "Elevated"
    assert roh.grade_for("brand-new").grade == "brand-new"
    assert roh.grade_for(None).grade == "—"
    assert roh.grade_text("none", show_relationship=True) == "None"
    assert "second-degree-like" in roh.grade_text(
        "second-degree-like", show_relationship=True
    )
    assert "degree" not in roh.grade_text("second-degree-like", show_relationship=False)


def test_viewer_badges_never_carry_the_relationship():
    viewer = roh.grade_badge("first-degree-like", show_relationship=False)
    curator = roh.grade_badge("first-degree-like", show_relationship=True)
    assert "first-degree" not in viewer["label"] + viewer["tooltip"]
    assert "first-degree-like" in curator["label"]


@pytest.mark.parametrize(
    "phi, degree",
    [(0.4, 0), (0.25, 1), (0.12, 2), (0.05, 3), (0.0, None), (-0.1, None)],
)
def test_kinship_degree(phi, degree):
    assert roh.kinship_degree(phi)[0] == degree
    assert roh.kinship_degree(None) == (None, "—")


def test_declared_relations():
    parents = {
        "A": (None, None),
        "B": (None, None),
        "C": ("A", "B"),
        "D": ("A", "B"),
        "E": ("A", "X"),
    }
    assert roh.declared_relation("A", "C", parents).label == "Parent–offspring"
    assert roh.declared_relation("C", "A", parents).label == "Parent–offspring"
    assert roh.declared_relation("C", "D", parents).label == "Full sibs"
    assert roh.declared_relation("C", "E", parents).label == "Half sibs"
    assert roh.declared_relation("A", "B", parents).parents_pair
    # Two founders share (None, None) parents; that is not a sibship.
    assert (
        roh.declared_relation("A", "X", {"A": (None, None), "X": (None, None)}) is None
    )


def test_kinship_rows_flag_only_declared_mismatches():
    kin = pl.DataFrame(
        {"IID1": ["A", "A", "C"], "IID2": ["B", "C", "Z"], "n_sites": [1, 1, 1],
         "ibs0": [0.05, 0.0, 0.0], "kinship": [0.07, 0.01, 0.25]}
    )  # fmt: skip
    rows = {r["pair"]: r for r in roh.kinship_rows(kin, {"C": ("A", "B")})}
    assert rows["A – B"]["parents_pair"] and not rows["A – B"]["mismatch"]
    assert rows["A – C"]["mismatch"]  # declared parent-offspring, unrelated
    assert not rows["C – Z"]["is_declared"] and not rows["C – Z"]["mismatch"]


def test_parents_from_froh_without_parent_columns():
    assert roh.parents_from_froh(pl.DataFrame({"IID": ["A"]})) == {"A": (None, None)}


def test_quantile_and_reference_values(app_config):
    assert roh.quantile([1, 2, 3, 4], 0.5) == 2.5
    assert roh.quantile([5, 1], 0) == 1 and roh.quantile([5, 1], 1) == 5
    assert roh.quantile([], 0.5) is None
    ref = roh.load_reference_froh("1.1.0")
    assert len(roh.reference_values(ref, "EUR", "froh_1.5mb")) == 30
    assert len(roh.reference_values(ref, roh.POOLED, "froh_1.5mb")) == 60
    bg = roh.reference_values(ref, "EUR", roh.F_BACKGROUND)
    assert bg[0] == pytest.approx(ref["f_short"][0] + ref["f_intermediate"][0])
    assert roh.reference_values(ref, "EUR", "not_a_column") == []
    assert roh.load_reference_hotspots("1.1.0")["chrom"].to_list() == ["1", "2"]


def test_reference_unconfigured_is_none():
    assert roh.load_reference_froh("1.1.0") is None
    assert roh.get_reference_roh_dir("1.1.0") is None


@pytest.mark.parametrize(
    "value, text",
    [
        (1, "1st"),
        (2, "2nd"),
        (3, "3rd"),
        (11, "11th"),
        (12, "12th"),
        (21, "21st"),
        (96.48, "96th"),
        (None, "—"),
    ],
)
def test_ordinal(value, text):
    assert roh.ordinal(value) == text


def test_flags():
    row = {
        "froh_category": "none",
        "homozygosity_pattern": "typical",
        "upd_candidate_chrom": None,
        "qc_pass": True,
    }
    assert not roh.is_flagged(row) and roh.flag_badges(row) == []
    upd = {**row, "upd_candidate_chrom": "8"}
    assert roh.is_flagged(upd) and roh.flag_badges(upd)[0]["label"] == "UPD? chr8"
    assert roh.is_flagged({**row, "qc_pass": False, "call_rate": 0.85})
    assert roh.is_flagged({**row, "homozygosity_pattern": "population"})


def test_to_float_and_metric_value():
    assert (
        roh.to_float("nan") is None
        and roh.to_float(True) is None
        and roh.to_float("x") is None
    )
    assert roh.metric_value(
        {"f_short": 0.1, "f_intermediate": 0.2}, roh.F_BACKGROUND
    ) == pytest.approx(0.3)
    assert (
        roh.metric_value({"f_short": None, "f_intermediate": 0.2}, roh.F_BACKGROUND)
        is None
    )
