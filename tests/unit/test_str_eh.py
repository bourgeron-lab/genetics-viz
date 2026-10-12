"""Short tandem repeats: discovery, reading and shaping (utils/str_eh.py)."""

import os

import polars as pl

from genetics_viz.utils import str_eh, str_plots


class Cohort:
    name = "DEMO"

    def __init__(self, data_dir):
        self.path = data_dir / "cohorts" / "DEMO"


def _family(data_dir, fid="1001"):
    d = str_eh.get_family_eh_dir(data_dir, fid)
    return str_eh.read_eh_table(str_eh.eh_path(d, fid, "jmsa50", "genotypes"))


def test_labels_newest_first_and_dotted(tmp_path):
    for i, label in enumerate(["old", "v1.2.cat"]):
        f = tmp_path / f"F.1.{label}.genotypes.tsv"
        f.write_text("x\n")
        os.utime(f, (1000 + i, 1000 + i))
    (tmp_path / "F.1.old.loci.tsv").write_text("x\n")
    (tmp_path / "OTHER.x.genotypes.tsv").write_text("x\n")
    assert str_eh.find_eh_labels(tmp_path, "F.1") == ["v1.2.cat", "old"]
    assert str_eh.find_eh_labels(tmp_path / "missing", "F.1") == []


def test_probes(data_dir):
    assert str_eh.probe_eh_data(data_dir, "1001")
    assert str_eh.probe_eh_data(data_dir, "AU-X-1")
    assert not str_eh.probe_eh_data(data_dir, "1002")
    assert str_eh.probe_cohort_eh_data(Cohort(data_dir))


def test_read_keeps_ids_as_text_and_dots_as_null(tmp_path):
    header = "family_id\tsample\tgenotype\tallele1\tallele2\tstatus"
    rows = [f"{100 + i}\tS{i}\t0/1\t10\t12\tnormal" for i in range(150)]
    rows += ["AU-ANG-FER-MPX3\tSx\t1\t9\t.\t."]
    path = tmp_path / "x.genotypes.tsv"
    path.write_text("\n".join([header, *rows]) + "\n")
    df = str_eh.read_eh_table(path)
    assert df.schema["family_id"] == pl.String
    assert df.schema["genotype"] == pl.String
    assert df.schema["allele2"] == pl.Int64
    assert df["allele2"][-1] is None and df["status"][-1] is None
    assert str_eh.read_eh_table(None) is None
    assert str_eh.read_eh_table(tmp_path / "nope.tsv") is None


def test_family_rows_without_cohort(data_dir):
    fam = _family(data_dir)
    members = str_eh.family_members(fam)
    assert [m.sample for m in members] == ["S1", "S2", "S3"]
    assert members[2].column == "S3 ♂●"

    rows = str_eh.family_locus_rows(fam, members, ["S1", "S2", "S3"])
    assert [r["locus"] for r in rows][:2] == ["HTT", "FXN"]  # most severe first
    htt = rows[0]
    assert htt["worst_key"] == "pathogenic"
    assert [b["label"] for b in htt["m2_badges"]] == ["18", "40·IR"]
    assert htt["m2_max"] == 40
    assert any(b["label"] == "↑ S3" for b in htt["flags_badges"])
    assert htt["founders_max"] is None  # no cohort context

    fxn = rows[1]
    assert fxn["worst_key"] == "carrier"
    assert not any(b["label"].startswith("↑") for b in fxn["flags_badges"])

    by_locus = {r["locus"]: r for r in rows}
    assert [b["label"] for b in by_locus["FMR1"]["m0_badges"]] == ["12"]  # male chrX
    assert by_locus["AFF3"]["worst_key"] == "unclassified"
    assert by_locus["AFF3"]["thresholds"] == "no thresholds"
    assert by_locus["HTT"]["thresholds"] == "N ≤26 · I ≥27 · P ≥36"


def test_family_rows_follow_the_selection(data_dir):
    fam = _family(data_dir)
    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(fam, members, ["S1", "S2"])
    htt = next(r for r in rows if r["locus"] == "HTT")
    assert htt["worst_key"] == "normal"
    assert "m2_badges" not in htt  # column ids stay stable: S3 is m2
    assert rows[0]["locus"] == "FXN"


def test_low_depth_is_shown(data_dir):
    fam = _family(data_dir, "AU-X-1")
    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(fam, members, ["S5", "S6"])
    arx = next(r for r in rows if r["locus"] == "ARX_2")
    assert "LowDepth" in [b["label"] for b in arx["m1_badges"]]


def test_cohort_context(data_dir):
    context = str_eh.load_cohort_context(Cohort(data_dir), "jmsa50")
    assert context.available and context.n_founders == 3
    assert context.genotypes.schema["family_id"] == pl.String

    fam = _family(data_dir)
    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(fam, members, ["S1", "S2", "S3"], context)
    fxn = next(r for r in rows if r["locus"] == "FXN")
    assert fxn["founders_max"] == 120 and fxn["cohort_flagged"] == 1
    assert "cohort pct" in fxn["m1_badges"][1]["tooltip"]

    assert str_eh.locus_distribution(context.genotypes, "FXN") == {12: 2, 120: 1}

    details = str_eh.member_details(fam, "HTT", context)
    assert {d["sample"] for d in details} == {"S1", "S2", "S3"}
    assert all(d["cohort_percentile"] is not None for d in details)


def test_missing_cohort_tables_degrade(data_dir):
    context = str_eh.load_cohort_context(Cohort(data_dir), "other_catalog")
    assert not context.available and context.n_founders is None


def test_cohort_table_cache_follows_the_file(data_dir):
    path = str_eh.eh_path(
        str_eh.get_cohort_eh_dir(Cohort(data_dir)), "DEMO", "jmsa50", "loci"
    )
    first = str_eh.read_eh_table_cached(path)
    assert str_eh.read_eh_table_cached(path) is first
    path.write_text(path.read_text() + path.read_text().splitlines()[1] + "\n")
    os.utime(path, (1, 1))
    assert len(str_eh.read_eh_table_cached(path)) == len(first) + 1


def test_common_in_cohort():
    row = {"n_samples": 500, "n_pathogenic": 96, "n_carrier": 0, "n_intermediate": 0}
    assert str_eh.common_in_cohort(row) == 96 / 500
    assert str_eh.common_in_cohort({**row, "n_pathogenic": 2}) is None
    assert str_eh.common_in_cohort(None) is None


def test_reviewer_path(data_dir):
    path = str_eh.reviewer_path(data_dir, "S3", "jmsa50", "HTT")
    assert path.exists() and path.name == "S3.jmsa50.HTT.svg"


def test_allele_class_matches_the_pipeline_order():
    thr = {"normal_max": 26, "intermediate_min": 27, "pathogenic_min": 36}
    assert [str_plots.allele_class(c, thr) for c in (20, 27, 36)] == [
        "normal",
        "intermediate",
        "pathogenic",
    ]
    gap = {"normal_max": 30, "intermediate_min": None, "pathogenic_min": 70}
    assert str_plots.allele_class(50, gap) == "above_normal"
    assert str_plots.allele_class(5, {}) == "unclassified"


def test_distribution_clips_a_far_threshold():
    option = str_plots.distribution_option(
        {10: 50, 14: 3},
        {"normal_max": None, "intermediate_min": None, "pathogenic_min": 1000},
        [("S1", [10, 16], "normal")],
    )
    assert option["xAxis"][1]["max"] < 100
    assert "pathogenic ≥ 1000" in option["title"]["subtext"]
    assert option["series"][1]["markLine"]["data"] == []
    assert len(option["series"][0]["data"]) == 2
