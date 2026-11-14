"""Genome QC (ghfc-ngs genomeqc step): discovery, reading, wording."""

import datagen
import polars as pl
import pytest

from genetics_viz.utils import genome_qc as gqc
from genetics_viz.utils import sample_qc as qc


class Cohort:
    def __init__(self, data_dir, name="DEMO"):
        self.name = name
        self.path = data_dir / "cohorts" / name


def test_no_tables_no_genome_qc(data_dir):
    # Default test data has sampleqc tables only: never read as genome QC.
    cohort = Cohort(data_dir)
    assert gqc.find_cohort_gqc_labels(gqc.get_cohort_gqc_dir(cohort), "DEMO") == []
    assert gqc.load_cohort_gqc(cohort) is None
    assert gqc.load_family_gqc(data_dir, "1001", cohort) is None
    assert not gqc.probe_cohort_gqc(cohort)
    assert not gqc.probe_family_gqc(data_dir, "1001", cohort)


def test_a_sampleqc_label_without_snapshot_is_not_genome_qc(tmp_path):
    # Strict: genomeqc always writes inputs.tsv; sampleqc's cohort tables share
    # the samples kind and an older run may have no snapshot.
    (tmp_path / "C.somalier_hg38.samples.tsv").write_text("sample\n")
    assert gqc.find_cohort_gqc_labels(tmp_path, "C") == []


def test_cohort_tables(data_dir):
    datagen.add_genomeqc_tables(data_dir)
    cohort = Cohort(data_dir)
    data = gqc.load_cohort_gqc(cohort)
    assert data.label == datagen.GQC_LABEL and data.source == "cohort"
    assert data.flagged
    assert set(data.samples.columns) >= set(gqc.COHORT_COLUMNS)
    assert data.samples["family_id"].dtype == pl.Utf8
    s3 = data.row_for("S3")
    assert s3["seq_class"] == "rl150" and s3["ebv_copies_per_cell"] == 152.3
    assert data.row_for("S1")["seq_class"] == "rl150_pcr"
    assert data.row_for("S6")["seq_class"] == "rl100"
    assert data.cell_lines() == {"S3"}
    # Sample QC's flag rows: one sample, level genome, low read as info,
    # most severe first.
    codes = [(f.code, f.severity, f.level) for f in data.flags]
    assert codes == [
        ("CELL_LINE", "high", "genome"),
        ("LOW_DEPTH", "medium", "genome"),
        ("EBV_TRACE", "info", "genome"),
    ]
    cell = data.flags_for("S3")[0]
    assert cell.keys == [("1001", "S3")]
    assert cell.ack_key == ("CELL_LINE", "1001", "S3", "")
    assert cell.info.category == "dna_source"
    assert cell.info.reader_label == "DNA from a cell line"
    # Sample QC still reads its own tables beside these.
    assert qc.load_cohort_qc(cohort).label == datagen.QC_LABEL


def test_family_from_the_cohort_tables(data_dir):
    datagen.add_genomeqc_tables(data_dir)
    cohort = Cohort(data_dir)
    fam = gqc.load_family_gqc(data_dir, "1001", cohort)
    assert fam.source == "cohort" and fam.prefix == "1001"
    assert sorted(fam.samples["sample"].to_list()) == ["S1", "S2", "S3"]
    assert [f.code for f in fam.flags] == ["CELL_LINE"]
    other = gqc.load_family_gqc(data_dir, "AU-X-1", cohort)
    assert {f.code for f in other.flags} == {"EBV_TRACE", "LOW_DEPTH"}
    assert gqc.probe_family_gqc(data_dir, "1001", cohort)


def test_family_without_cohort_tables(data_dir):
    # A standalone family, or a cohort not gathered yet: the family's
    # variants joined to its members' sample tables, no class, no flag.
    datagen.add_genomeqc_tables(data_dir, cohort=None)
    assert gqc.load_cohort_gqc(Cohort(data_dir)) is None
    fam = gqc.load_family_gqc(data_dir, "1001", Cohort(data_dir))
    assert fam.source == "family" and not fam.flagged and fam.flags == []
    s3 = fam.row_for("S3")
    assert s3["ebv_copies_per_cell"] == 152.3 and s3["het_hom_ratio"] == 1.5293
    assert s3["seq_class"] is None and s3["flags"] is None
    assert set(fam.samples.columns) >= set(gqc.COHORT_COLUMNS)
    assert gqc.probe_family_gqc(data_dir, "1001")
    assert gqc.load_family_gqc(data_dir, "1002") is None


def test_sample_tables_and_stats(data_dir):
    datagen.add_genomeqc_tables(data_dir)
    directory = gqc.get_sample_gqc_dir(data_dir, "S3")
    assert gqc.find_sample_gqc_labels(directory, "S3") == [datagen.GQC_LABEL]
    path = gqc.sample_stats_path(data_dir, "S3")
    assert path is not None and path.name == "S3.gqc_v1.stats.txt"
    assert gqc.sample_stats_path(data_dir, "S4") is None


def test_read_keeps_ids_as_text_and_adds_missing_columns(tmp_path):
    rows = [f"S{i}\t{88171200000 + i}\t." for i in range(50)] + ["Sx\tAU-X\tCELL_LINE"]
    path = tmp_path / "x.samples.tsv"
    path.write_text("\n".join(["sample\tfamily_id\tflags", *rows]) + "\n")
    df = gqc.read_gqc_table(path, gqc.COHORT_COLUMNS)
    assert df["family_id"][0] == "88171200000"
    assert df["flags"][0] is None
    assert "ebv_copies_per_cell" in df.columns


def test_metric_wording():
    ebv = gqc.metric_info("ebv_copies_per_cell")
    assert ebv.log and ebv.lines == [0.01, 1.0] and ebv.group == "dna_source"
    assert ebv.fmt(0.0174) == "0.017"
    assert gqc.metric_info("duplicate_share").fmt(0.1202) == "12.0%"
    assert gqc.metric_info("n_indel").fmt(587299) == "587,299"
    assert gqc.metric_info("n_indel").side == "both"
    assert gqc.metric_info("chr20_depth").lines == [20.0]
    assert gqc.metric_info("private_snvs").plot is False
    assert gqc.metric_info("library").fmt("pcr") == "pcr"
    assert gqc.metric_info("ti_tv").fmt(None) == "—"
    unknown = gqc.metric_info("new_metric")
    assert unknown.label == "new_metric" and unknown.fmt(1.5) == "1.50"
    assert gqc.seq_class_info("rl150_pcr")["label"] == "150 bp, PCR"
    assert gqc.seq_class_info(None)["key"] == "unknown"
    assert gqc.seq_class_info("rl75")["label"] == "rl75"


@pytest.mark.parametrize("group", ["dna_source", "alignment", "variants"])
def test_every_metric_has_a_group(group):
    assert gqc.metrics(group)
    assert {m.group for m in gqc.metrics()} <= set(gqc.GROUPS)


def test_family_without_sample_tables(data_dir):
    # Members with no CRAM have no sample table: the metrics are null.
    datagen.add_genomeqc_tables(data_dir, cohort=None)
    for sid in ("S1", "S2", "S3"):
        (gqc.get_sample_gqc_dir(data_dir, sid) / f"{sid}.gqc_v1.tsv").unlink()
    fam = gqc.load_family_gqc(data_dir, "1001")
    assert fam.row_for("S3")["chr20_depth"] is None
    assert fam.row_for("S3")["n_indel"] == 587299
