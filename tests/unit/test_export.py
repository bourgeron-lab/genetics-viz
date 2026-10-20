"""TSV exports of the cohort and family tabs (utils/export.py and the builders)."""

import io
import zipfile
from datetime import date

import datagen
import pytest

from genetics_viz.utils import export
from genetics_viz.utils import sample_qc as qc
from genetics_viz.utils.ancestry import (
    ancestry_export,
    find_ancestry_files,
    find_cohort_ancestry_files,
    pgs_export,
)
from genetics_viz.utils.roh import (
    find_roh_files,
    froh_export,
    kinship_export,
    read_roh_table,
    segments_export,
)
from genetics_viz.utils.sample_qc_acks import load_acks
from genetics_viz.utils.sample_qc_export import (
    flags_export,
    pairs_export,
    samples_export,
)
from genetics_viz.utils.str_eh import (
    FLAGGED,
    eh_path,
    genotypes_export,
    get_family_eh_dir,
    load_cohort_context,
    read_eh_table,
)


class Cohort:
    def __init__(self, data_dir):
        self.name = "DEMO"
        self.path = data_dir / "cohorts" / "DEMO"


def _people(*rows):
    """``export.person`` dicts from ``(FID, IID, sex, phenotype)`` tuples."""
    return [export.person(*r) for r in rows]


def _demo_people(*iids):
    return _people(
        *[(m[0], m[1], m[4], m[5]) for m in datagen.DEMO_PEDIGREE if m[1] in iids]
    )


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------


def test_rows_to_tsv_writes_missing_as_empty_and_joins_lists():
    text = export.rows_to_tsv(
        [{"a": None, "b": ["x", "y"], "c": float("nan"), "d": 1.5}],
        ["a", "b", "c", "d"],
        headers=["A", "B", "C", "D"],
    )
    assert text == "A\tB\tC\tD\n\tx, y\t\t1.5\n"


def test_is_export_column_drops_display_state():
    assert export.is_export_column("VEP_SYMBOL")
    for key in ("actions", "GeneBadges", "flags_badges", "Diagnostic_badge"):
        assert not export.is_export_column(key)


def test_export_filename_is_dated_and_safe():
    name = export.export_filename("AU X/1", "roh", "genes", tag="apgs_b1.1.0_dp6gq15")
    assert name == f"AU_X_1.roh.genes.apgs_b1.1.0_dp6gq15.{date.today()}.tsv"


def test_zip_texts_round_trip():
    data = export.zip_texts({"a.tsv": "x\n1\n", "b.tsv": "y\n2\n"})
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        assert archive.namelist() == ["a.tsv", "b.tsv"]
        assert archive.read("b.tsv").decode() == "y\n2\n"


def test_people_are_normalised_and_counted_once():
    individuals = [
        {"FID": "1002", "Sample ID": "S4", "Sex": "2.0", "Phenotype": "-9.0"},
        {"FID": "1002", "Sample ID": "S4", "Sex": "2.0", "Phenotype": "-9.0"},
        {"FID": "1001", "Sample ID": "S1", "Sex": "-", "Phenotype": "2"},
    ]
    assert export.people_from_individuals(individuals) == [
        {"FID": "1002", "IID": "S4", "sex": "2", "phenotype": "-9"},
        {"FID": "1001", "IID": "S1", "sex": "-9", "phenotype": "2"},
    ]


# ---------------------------------------------------------------------------
# Ancestry and PGS
# ---------------------------------------------------------------------------


def test_cohort_ancestry_export_has_labels_pcs_and_admixture(data_dir):
    files = find_cohort_ancestry_files(Cohort(data_dir))
    df = ancestry_export(files, _demo_people("S1", "S5"))
    assert df["IID"].to_list() == ["S1", "S5"]
    assert df["FID"].to_list() == ["1001", "AU-X-1"]
    assert df.columns[:6] == ["FID", "IID", "sex", "phenotype", "region", "population"]
    assert {"PC1", "PC20", "snp_coverage", "admix_EUR", "admix_CSA"} <= set(df.columns)
    assert "family_id" not in df.columns
    assert df["admix_EUR"].to_list() == [0.8, 0.8]


def test_family_ancestry_export_keeps_only_people_with_pcs(data_dir):
    files = find_ancestry_files(data_dir, "1001")
    people = _demo_people("S1", "S3") + _people(("1001", "NOT_SEQUENCED", "1", "1"))
    df = ancestry_export(files, people)
    assert df["IID"].to_list() == ["S1", "S3"]
    assert df["sex"].to_list() == ["1", "1"] and df["phenotype"].to_list() == ["1", "2"]


def test_pgs_export_follows_the_region_selection(data_dir):
    files = find_cohort_ancestry_files(Cohort(data_dir))
    people = _demo_people("S1", "S2", "S4")
    df = pgs_export(files, people)
    assert df.columns == ["FID", "IID", "sex", "phenotype", "region", *datagen.TRAITS]
    assert df["region"].to_list() == ["EUR"] * 3
    assert pgs_export(files, people, regions={"CSA"}).height == 0
    assert pgs_export(files, people, regions={"EUR"}).height == 3


# ---------------------------------------------------------------------------
# ROH
# ---------------------------------------------------------------------------


def _froh_rows(data_dir, fid="1001"):
    files = find_roh_files(data_dir, fid)
    return files, read_roh_table(files.path("froh"))


@pytest.mark.parametrize("curator", [True, False])
def test_froh_export_wording_follows_the_role(data_dir, curator):
    _, froh = _froh_rows(data_dir)
    people = {p["IID"]: p for p in _demo_people("S1", "S2", "S3")}
    paired = [(r, people[r["IID"]]) for r in froh.to_dicts()]
    rows, columns = froh_export(paired, show_relationship=curator)
    assert {"froh_1.5mb", "pct_froh_1.5mb_vs_ref_group", "f_long", "pattern",
            "grade", "flags", "f_short+intermediate"} <= set(columns)  # fmt: skip
    hidden = {"froh_category", "parents_kinship", "expected_f"}
    assert (hidden <= set(columns)) is curator
    assert not (hidden & set(columns)) or curator
    child = next(r for r in rows if r["IID"] == "S3")
    text = export.rows_to_tsv([child], columns)
    assert ("first-cousin" in text or "degree" in text) is curator


def test_family_segments_and_kinship_exports(data_dir):
    files, froh = _froh_rows(data_dir)
    segments = segments_export(read_roh_table(files.path("roh")), "1001", ["S3"])
    assert segments.columns[:2] == ["FID", "IID"]
    assert segments["IID"].to_list() == ["S3"]

    kinship = read_roh_table(files.path("kinship"))
    members = ["S1", "S2", "S3"]
    curator_rows, columns = kinship_export(kinship, froh, members, True)
    reader_rows, _ = kinship_export(kinship, froh, members, False)
    assert "observed" in columns
    # The parents' pair is shown to curators only.
    assert {(r["IID1"], r["IID2"]) for r in curator_rows} == {
        ("S1", "S2"), ("S1", "S3"), ("S2", "S3"),
    }  # fmt: skip
    assert {(r["IID1"], r["IID2"]) for r in reader_rows} == {("S1", "S3"), ("S2", "S3")}


# ---------------------------------------------------------------------------
# STR
# ---------------------------------------------------------------------------


def test_cohort_flagged_genotypes_export(data_dir):
    context = load_cohort_context(Cohort(data_dir), datagen.EH_LABEL)
    df = genotypes_export(
        context.genotypes, _demo_people("S2", "S3", "S1"), statuses=FLAGGED
    )
    assert df.columns[:4] == ["FID", "IID", "sex", "phenotype"]
    assert {"locus", "allele1", "allele2", "status", "cohort_percentile"} <= set(
        df.columns
    )
    assert set(df["status"].to_list()) <= set(FLAGGED)
    assert ("S3", "HTT") in set(zip(df["IID"].to_list(), df["locus"].to_list()))
    assert "S1" not in df["IID"].to_list()
    # Pedigree codes, not the pipeline's words.
    assert set(df["sex"].to_list()) <= {"1", "2", "-9"}


def test_family_genotypes_export_adds_cohort_and_curation(data_dir):
    fam = read_eh_table(
        eh_path(get_family_eh_dir(data_dir, "1001"), "1001", datagen.EH_LABEL,
                "genotypes")
    )  # fmt: skip
    context = load_cohort_context(Cohort(data_dir), datagen.EH_LABEL)
    df = genotypes_export(
        fam,
        _demo_people("S2", "S3"),
        loci={"HTT", "FXN"},
        context=context,
        validations={("HTT", "S3"): "present"},
    )
    assert set(df["locus"].to_list()) == {"HTT", "FXN"}
    assert set(df["IID"].to_list()) == {"S2", "S3"}
    assert df["cohort_percentile"].null_count() < df.height
    assert df.columns[-1] == "validation"
    by_key = dict(
        zip(zip(df["IID"].to_list(), df["locus"].to_list()), df["validation"])
    )
    assert by_key[("S3", "HTT")] == "present" and by_key[("S2", "HTT")] is None


# ---------------------------------------------------------------------------
# Sample QC
# ---------------------------------------------------------------------------


def _cohort_qc(data_dir):
    return qc.load_cohort_qc(Cohort(data_dir)), load_acks(data_dir)


@pytest.mark.parametrize("curator", [True, False])
def test_flags_export_carries_the_acknowledgements(data_dir, curator):
    data, acks = _cohort_qc(data_dir)
    rows, columns = flags_export(data.flags, acks, show_detail=curator)
    assert ("evidence" in columns) is curator
    assert columns[-4:] == ["ack_status", "ack_note", "ack_user", "ack_date"]
    by_code = {r["code"]: r for r in rows}
    cross = by_code["CRYPTIC_RELATED"]
    assert cross["ack_status"] == "explained" and cross["ack_date"] == "2026-09-03"
    # A relatedness note names a relationship; a sex note does not.
    assert (cross["ack_note"] == "(hidden)") is not curator
    assert by_code["SEX_MISMATCH"]["ack_note"].startswith("Pedigree entry error")
    # The withdrawn acknowledgement leaves the flag counting.
    assert by_code["POSSIBLE_CONTAMINATION"]["ack_status"] == ""
    # Acknowledged flags come after the others.
    assert [bool(r["ack_status"]) for r in rows] == sorted(
        bool(r["ack_status"]) for r in rows
    )
    if not curator:
        assert cross["issue"] == qc.flag_info("CRYPTIC_RELATED").reader_label


@pytest.mark.parametrize("curator", [True, False])
def test_samples_export_masks_karyotype_and_splits_acknowledged(data_dir, curator):
    data, acks = _cohort_qc(data_dir)
    people = {(p["FID"], p["IID"]): p for p in _demo_people("S3", "S5", "S6")}
    rows, columns = samples_export(
        data.sample_rows(), people, data.flags, acks, curator
    )
    by_iid = {r["IID"]: r for r in rows}
    assert set(by_iid) == {"S3", "S5", "S6"}
    assert by_iid["S5"]["karyotype"] == ("XXY" if curator else "Atypical")
    assert by_iid["S6"]["flags"] == ""
    assert "(explained)" in by_iid["S6"]["acknowledged_flags"]
    assert by_iid["S6"]["worst_severity"] == ""
    assert by_iid["S3"]["worst_severity"] == "high"
    assert {"depth_mean", "contamination_charr", "x_copies"} <= set(columns)


def test_pairs_export_hides_undeclared_pairs_from_readers(data_dir):
    data, acks = _cohort_qc(data_dir)
    curator_rows, curator_columns = pairs_export(
        data.pair_rows(), data.flags, acks, True, first_family="AU-X-1"
    )
    reader_rows, reader_columns = pairs_export(
        data.pair_rows(), data.flags, acks, False
    )
    assert "observed_relationship" in curator_columns
    assert "observed_relationship" not in reader_columns

    def pairs(rows):
        return {(r["sample_a"], r["sample_b"]) for r in rows}

    # The cross pair, AU-X-1's member first, and the undeclared parents' pair.
    assert {("S5", "S2"), ("S1", "S2")} <= pairs(curator_rows)
    assert pairs(reader_rows) == {("S1", "S3"), ("S2", "S3"), ("S5", "S6")}
    cross = next(
        r for r in curator_rows if r["sample_a"] == "S5" and r["sample_b"] == "S2"
    )
    assert cross["FID_a"] == "AU-X-1"
    assert "(explained)" in cross["acknowledged_flags"]
    mismatch = next(r for r in reader_rows if r["sample_a"] == "S1")
    assert mismatch["flags"] == qc.flag_info("NOT_PARENT").reader_label
