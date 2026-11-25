"""Aneuploidy validations (utils/aneuploidy.py), the flags they acknowledge,
and what the report and Diagnostics panel say about them."""

import datagen
import pytest

from genetics_viz.components.validation_loader import load_validation_map
from genetics_viz.utils import aneuploidy as an
from genetics_viz.utils import cohort_svs
from genetics_viz.utils import sample_qc as qc
from genetics_viz.utils import sample_qc_acks as acks
from genetics_viz.utils.diagnostic_labels import describe_aneuploidy
from genetics_viz.utils.validation_files import validation_path

#: AU-X-1: S5 male, XXY in Sample QC; S6 declared female, XY.
AU_MEMBERS = [("S5", "1"), ("S6", "2")]
#: 1001: S1 father, S2 mother, S3 son.
MEMBERS_1001 = [("S1", "1"), ("S2", "2"), ("S3", "1")]


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(an, "can_write", lambda: True)


def _save(data_dir, **kw):
    args = dict(
        family_id="AU-X-1",
        sample="S5",
        chrom="XY",
        user="curator",
        validation="present",
    )
    args.update(kw)
    an.save_aneuploidy_validation(data_dir, **args)


def _karyotype_flag(sample="S5", fid="AU-X-1"):
    return qc.FlagRow("family", fid, sample, None, "KARYOTYPE_XXY", "high", "")


def test_keys():
    assert an.aneuploidy_key("21") == an.aneuploidy_key("chr21") == "chr21:ANEUPLOIDY"
    # X and Y are one record: a karyotype such as XXY is about both.
    assert an.aneuploidy_key("chrX") == an.aneuploidy_key("Y") == "chrXY:ANEUPLOIDY"
    assert an.chrom_of_key("chrXY:ANEUPLOIDY") == "XY"
    assert an.chrom_label("XY") == "sex chromosomes"
    assert an.chrom_label("21") == "chr21"


@pytest.mark.parametrize(
    "chrom, direction, observed, declared, expected",
    [
        ("21", "gain", "XY", "", "47,XY,+21"),
        ("18", "loss", "", "XX", "45,XX,-18"),
        # An atypical observed karyotype is not used for an autosome.
        ("21", "gain", "XXY", "XY", "47,XY,+21"),
        ("XY", "", "XXY", "XY", "47,XXY"),
        ("XY", "", "X0", "XX", "45,X"),
        ("21", "", "XY", "", ""),
        ("21", "gain", "", "", ""),
    ],
)
def test_suggested_karyotype(chrom, direction, observed, declared, expected):
    assert an.suggested_karyotype(chrom, direction, observed, declared) == expected


def test_viewers_and_bad_values_are_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(an, "can_write", lambda: False)
    with pytest.raises(PermissionError):
        _save(tmp_path)
    monkeypatch.setattr(an, "can_write", lambda: True)
    for bad in (
        dict(validation="maybe"),
        dict(origin="grandmaternal"),
        dict(direction="sideways"),
        dict(mosaic="partly"),
        dict(chrom=""),
        dict(user=""),
    ):
        with pytest.raises(ValueError):
            _save(tmp_path, **bad)
    assert not an.aneuploidy_validation_file(tmp_path).exists()


def test_shared_loader_reads_the_file(tmp_path, curator):
    _save(tmp_path, origin="de novo", karyotype="47,XXY ", comment="array")
    path = an.aneuploidy_validation_file(tmp_path)
    assert path == validation_path(tmp_path, "Aneuploidy")
    header, row = path.read_text().splitlines()
    assert header.split("\t")[-3:] == ["Direction", "Mosaic", "Karyotype"]
    assert row.split("\t")[-1] == "47,XXY"
    vmap = load_validation_map(path, "AU-X-1")
    validation, origin, comment, *_ = vmap[("chrXY:ANEUPLOIDY", "S5")][0]
    assert (validation, origin, comment) == ("present", "de novo", "array")


def test_sex_chromosomes_from_sample_qc(data_dir):
    review = an.load_family_review(data_dir, "AU-X-1", AU_MEMBERS)
    by_sample = {c.sample: c for c in review.candidates}
    s5 = by_sample["S5"]
    assert (s5.chrom, s5.karyotype, s5.label) == ("XY", "XXY", "XXY")
    assert [f.code for f in s5.flags] == ["KARYOTYPE_XXY"]
    assert s5.direction == "gain" and s5.suggested_karyotype == "47,XXY"
    # S6 is XY but declared female: listed, as the report lists it.
    assert by_sample["S6"].karyotype == "XY"
    assert an.count_to_review(review.candidates) == 2


def test_autosomes_from_wisecondorx(data_dir):
    assert an.load_family_review(data_dir, "1001", MEMBERS_1001).candidates == []
    datagen.add_report_case(data_dir)
    review = an.load_family_review(data_dir, "1001", MEMBERS_1001)
    (c,) = review.candidates
    assert (c.sample, c.chrom, c.label, c.direction) == ("S3", "21", "+21", "gain")
    assert c.events[0].text == "+21 (100%)"
    assert c.suggested_karyotype == "47,XY,+21"
    assert [s["chr"] for s in review.segments("S3", "21")] == ["chr21"]


def test_a_validation_keeps_its_chromosome_listed(data_dir, curator):
    _save(data_dir, family_id="1001", sample="S2", chrom="18", validation="absent")
    review = an.load_family_review(data_dir, "1001", MEMBERS_1001)
    (c,) = review.candidates
    assert (c.sample, c.chrom, c.status, c.label) == ("S2", "18", "absent", "chr18")
    assert an.count_to_review(review.candidates) == 0


@pytest.mark.parametrize(
    "validation, status",
    [("present", "confirmed"), ("absent", "explained"), ("uncertain", None)],
)
def test_sex_chromosome_validation_acknowledges_the_flag(
    tmp_path, curator, validation, status
):
    _save(tmp_path, validation=validation, karyotype="47,XXY", comment="array")
    ack = acks.load_acks(tmp_path).get(_karyotype_flag())
    if status is None:
        assert ack is None
        return
    assert (ack.status, ack.source) == (status, "aneuploidy")
    assert ack.note == f"Sex chromosomes validated {validation} (47,XXY): array"


def test_only_the_sex_chromosomes_of_that_sample_acknowledge(tmp_path, curator):
    _save(tmp_path, chrom="21")
    assert _karyotype_flag() not in acks.load_acks(tmp_path)
    _save(tmp_path)
    standing = acks.load_acks(tmp_path)
    assert _karyotype_flag() in standing
    assert _karyotype_flag(sample="S6") not in standing
    assert _karyotype_flag(fid="1001") not in standing


def _aneuploidy_section(data_dir):
    from genetics_viz.models import DataStore
    from genetics_viz.utils.family_report import build_family_report

    store = DataStore(data_dir)
    store.load()
    cohort = store.cohorts["DEMO"]
    report = build_family_report(
        data_dir, "1001", cohort.families["1001"], cohort, "curator"
    )
    return next(s for s in report.sections if s.key == "aneuploidy")


def test_report_shows_the_validation(data_dir, app_config, curator):
    datagen.add_report_case(data_dir)
    section = _aneuploidy_section(data_dir)
    s3 = next(r for r in section.rows if r["iid"] == "S3")
    assert s3["validation"] == "+21: not reviewed"
    assert "S3: +21 (100%)" in section.glance

    _save(
        data_dir,
        family_id="1001",
        sample="S3",
        chrom="21",
        validation="absent",
        comment="artefact",
    )
    section = _aneuploidy_section(data_dir)
    s3 = next(r for r in section.rows if r["iid"] == "S3")
    assert s3["validation"] == "+21: absent"
    # Validated absent: no longer reported as found.
    assert "S3" not in section.glance


def test_diagnostic_wording_hides_the_karyotype_from_readers(tmp_path, curator):
    _save(tmp_path, karyotype="47,XXY")
    entry = {"Variant": "chrXY:ANEUPLOIDY", "Sample": "S5", "Impact": "Sample QC"}
    text, tooltip = describe_aneuploidy(tmp_path, "AU-X-1", entry)
    assert text == "47,XXY" and "S5: validated present, 47,XXY" in tooltip
    text, tooltip = describe_aneuploidy(tmp_path, "AU-X-1", entry, show_detail=False)
    assert text == "sex chromosomes" and "XXY" not in tooltip


def test_cohort_rows_word_each_event_validation():
    rows = [
        {
            "chr": "chr21",
            "start": 10_500_001,
            "end": 46_700_000,
            "ratio": 0.58,
            "FID": "1001",
            "sample": "S3",
        }
    ]
    (row,) = cohort_svs.aneuploidy_rows(rows, 0.75)
    assert row["validation"] == "+21: not reviewed"
    (row,) = cohort_svs.aneuploidy_rows(rows, 0.75, {("1001", "S3", "21"): "present"})
    assert row["validation"] == "+21: present"
