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
    # Most severe first: high, medium, info, then by status.
    assert [r["locus"] for r in rows][:3] == ["HTT", "FMR1", "FXN"]
    htt = rows[0]
    assert htt["worst_key"] == "pathogenic"
    assert [b["label"] for b in htt["severity_badges"]] == ["High"]
    assert [b["label"] for b in htt["m2_badges"]] == ["18", "40"]
    assert htt["m2_max"] == 40
    assert any(b["label"] == "↑ S3" for b in htt["flags_badges"])
    assert htt["founders_max"] is None  # no cohort context

    fmr1 = rows[1]
    assert fmr1["worst_key"] == "premutation"
    assert [b["label"] for b in fmr1["severity_badges"]] == ["Medium"]
    assert [b["label"] for b in fmr1["m0_badges"]] == ["12"]  # male chrX
    assert fmr1["thresholds"] == "N ≤40 · I ≥41 · PM ≥55 · P ≥201"

    fxn = rows[2]
    assert fxn["worst_key"] == "carrier"
    assert [b["label"] for b in fxn["severity_badges"]] == ["Info"]
    assert not any(b["label"].startswith("↑") for b in fxn["flags_badges"])

    by_locus = {r["locus"]: r for r in rows}
    assert by_locus["AFF3"]["worst_key"] == "unclassified"
    assert by_locus["AFF3"]["severity_badges"] == []
    assert by_locus["AFF3"]["thresholds"] == "no thresholds"
    assert by_locus["HTT"]["thresholds"] == "N ≤26 · I ≥27 · P ≥36"
    assert by_locus["TCF4"]["thresholds"] == "N ≤40 · P ≥50 · risk"


def test_family_rows_follow_the_selection(data_dir):
    fam = _family(data_dir)
    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(fam, members, ["S1", "S2"])
    htt = next(r for r in rows if r["locus"] == "HTT")
    assert htt["worst_key"] == "normal"
    assert "m2_badges" not in htt  # column ids stay stable: S3 is m2
    assert rows[0]["locus"] == "FMR1"


def test_low_depth_is_shown(data_dir):
    fam = _family(data_dir, "AU-X-1")
    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(fam, members, ["S5", "S6"])
    arx = next(r for r in rows if r["locus"] == "ARX_2")
    assert "LowDepth" in [b["label"] for b in arx["m1_badges"]]


def test_low_confidence_lowers_the_severity_and_marks_the_allele(data_dir):
    """S6's ARX_2 is pathogenic, but no read spans its GCN tract: medium."""
    fam = _family(data_dir, "AU-X-1")
    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(fam, members, ["S5", "S6"])
    arx = next(r for r in rows if r["locus"] == "ARX_2")
    assert arx["worst_key"] == "pathogenic"
    assert [b["label"] for b in arx["severity_badges"]] == ["Medium"]
    allele = arx["m1_badges"][1]
    assert allele["label"] == "22·F?"
    assert "low confidence" in allele["tooltip"]
    # A risk allele is never more than info, whatever its status.
    tcf4 = next(r for r in rows if r["locus"] == "TCF4")
    assert tcf4["worst_key"] == "pathogenic"
    assert [b["label"] for b in tcf4["severity_badges"]] == ["Info"]


def test_cohort_context(data_dir):
    context = str_eh.load_cohort_context(Cohort(data_dir), "jmsa50")
    assert context.available and context.n_founders == 3
    assert context.genotypes.schema["family_id"] == pl.String

    fam = _family(data_dir)
    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(fam, members, ["S1", "S2", "S3"], context)
    fxn = next(r for r in rows if r["locus"] == "FXN")
    assert fxn["founders_max"] == 120 and fxn["cohort_flagged"] == 0  # info
    fmr1 = next(r for r in rows if r["locus"] == "FMR1")
    assert fmr1["cohort_flagged"] == 1
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
    row = {"n_samples": 500, "n_pathogenic": 96, "n_high": 10, "n_medium": 16,
           "n_info": 70}  # fmt: skip
    assert str_eh.common_in_cohort(row) == 26 / 500
    assert str_eh.common_in_cohort({**row, "n_medium": 2}) is None
    assert str_eh.common_in_cohort(None) is None
    # A loci table written before severities counts the flagged statuses.
    legacy = {"n_samples": 500, "n_pathogenic": 96, "n_carrier": 0,
              "n_intermediate": 0}  # fmt: skip
    assert str_eh.common_in_cohort(legacy) == 96 / 500


def test_reviewer_path(data_dir):
    path = str_eh.reviewer_path(data_dir, "S3", "jmsa50", "HTT")
    assert path is not None and path.name == "S3.jmsa50.HTT.svg"
    assert path.parent.name == "reviewer"
    assert str_eh.reviewer_path(data_dir, "S3", "jmsa50", "FXN") is None
    assert str_eh.reviewer_reviewed(data_dir, "S3", "jmsa50")
    assert not str_eh.reviewer_reviewed(data_dir, "S3", "other")


def test_reviewer_path_maps_a_variant_to_its_locus(data_dir):
    """REViewer names images by LocusId: FXN_A and FXN_B are drawn as FXN."""
    rdir = str_eh.get_reviewer_dir(data_dir, "S3")
    for locus in ("FXN", "F"):
        (rdir / f"S3.jmsa50.{locus}.svg").write_text("<svg/>")
    path = str_eh.reviewer_path(data_dir, "S3", "jmsa50", "FXN_A")
    assert path is not None and path.name == "S3.jmsa50.FXN.svg"
    assert str_eh.reviewer_path(data_dir, "S3", "jmsa50", "FXNX") is None
    assert str_eh.reviewer_path(data_dir, "S3", "other", "FXN_A") is None


def test_allele_class_matches_the_pipeline_order():
    thr = {"normal_max": 26, "intermediate_min": 27, "pathogenic_min": 36}
    assert [str_plots.allele_class(c, thr) for c in (20, 27, 36)] == [
        "normal",
        "intermediate",
        "pathogenic",
    ]
    fmr1 = {"normal_max": 40, "intermediate_min": 41, "premutation_min": 55,
            "pathogenic_min": 201}  # fmt: skip
    assert [str_plots.allele_class(c, fmr1) for c in (45, 60, 250)] == [
        "intermediate",
        "premutation",
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


def test_restrict_to_keys_by_family_and_sample(data_dir):
    context = str_eh.load_cohort_context(Cohort(data_dir), "jmsa50")
    sub = str_eh.restrict_to(context.genotypes, [("1001", "S3"), ("9999", "S1")])
    assert set(sub["sample"].unique()) == {"S3"}
    assert len(str_eh.restrict_to(context.genotypes, [])) == 0


def test_cohort_rows_follow_the_selection(data_dir):
    context = str_eh.load_cohort_context(Cohort(data_dir), "jmsa50")
    everyone = [("1001", "S1"), ("1001", "S2"), ("1001", "S3"), ("AU-X-1", "S5"),
                ("AU-X-1", "S6")]  # fmt: skip
    loci = str_eh.cohort_locus_rows(context, everyone)
    assert [r["locus"] for r in loci[:3]] == ["HTT", "ARX_2", "FMR1"]
    htt, arx, fmr1 = loci[:3]
    assert (htt["n_high"], htt["n_pathogenic"], htt["founders_max"]) == (1, 1, 12)
    assert (arx["n_high"], arx["n_medium"], arx["n_pathogenic"]) == (0, 1, 1)
    assert (fmr1["n_medium"], fmr1["n_premutation"]) == (1, 1)
    fxn = next(r for r in loci if r["locus"] == "FXN")
    assert (fxn["n_carrier"], fxn["n_flagged"]) == (1, 0)

    flagged = str_eh.cohort_flagged_rows(context, everyone)
    assert [(r["sample"], r["locus"], r["status"], r["severity"]) for r in flagged] == [
        ("S3", "HTT", "Pathogenic", "High"),
        ("S6", "ARX_2", "Pathogenic", "Medium"),
        ("S2", "FMR1", "Premutation", "Medium"),
    ]
    assert flagged[0]["allele_values"] == [18, 40]

    parents = [k for k in everyone if k[1] != "S3"]
    assert str_eh.cohort_locus_rows(context, parents)[0]["locus"] == "ARX_2"
    assert [r["sample"] for r in str_eh.cohort_flagged_rows(context, parents)] == [
        "S6",
        "S2",
    ]


def test_tables_without_severity_flag_by_status(data_dir):
    """A family not yet rebuilt by the pipeline keeps the old status rule."""
    context = str_eh.load_cohort_context(Cohort(data_dir), "jmsa50")
    legacy = str_eh.CohortContext(
        context.name,
        context.label,
        context.loci.drop("n_high", "n_medium", "n_info"),
        context.genotypes.drop("severity"),
    )
    everyone = [("1001", "S2"), ("1001", "S3")]
    flagged = str_eh.cohort_flagged_rows(legacy, everyone)
    assert [(r["sample"], r["locus"]) for r in flagged] == [
        ("S3", "HTT"),
        ("S2", "FXN"),
    ]
    assert flagged[0]["severity_badges"] == []
    loci = {r["locus"]: r for r in str_eh.cohort_locus_rows(legacy, everyone)}
    assert loci["FXN"]["n_flagged"] == 1 and loci["FXN"]["n_high"] is None
    assert loci["FMR1"]["n_flagged"] == 0  # premutation was not a flagged status
    assert str_eh.is_flagged({"status": "carrier"})
    assert not str_eh.is_flagged({"status": "carrier", "severity": "info"})
