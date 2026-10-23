"""The cross-cohort gene index and the gene/variant page rows built from it."""

import os
from pathlib import Path

import polars as pl
import pytest

import datagen
from genetics_viz.config_model import DataDirectoryConfig
from genetics_viz.utils import data, gene_index, gene_rows
from genetics_viz.utils.gene_scoring import GeneScorer
from genetics_viz.utils.genetrek import GenetrekSnapshot
from genetics_viz.utils.user_preferences import GenePrefs
from genetics_viz.utils.sample_counts import read_cohort_wombat

SHANK3_DE_NOVO = "chr22:50674000:G:A"


@pytest.fixture
def indexed(app_config, data_dir: Path) -> Path:
    """The synthetic data directory, plus cohort OTHER, registered and indexed."""
    datagen.add_gene_case(data_dir)
    data.init_all_data_stores([DataDirectoryConfig(path=str(data_dir))])
    gene_index.refresh_all(wait=True)
    return data_dir


def _manifest(data_dir: Path) -> dict:
    return gene_index._manifests[str(data_dir)]


def _scorer() -> GeneScorer:
    return GeneScorer(GenetrekSnapshot(), GenePrefs(["#c62828"], {}))


# ---------------------------------------------------------------------------
# Reading the sources
# ---------------------------------------------------------------------------


def test_read_cohort_wombat_reads_only_the_named_columns(data_dir: Path):
    path = data_dir / "cohorts/DEMO/wombat/DEMO.rare.vep1.default.results.tsv"
    df = read_cohort_wombat(path, columns=["sample", "#CHROM", "absent"])
    # File order, whatever the order asked; a column the table lacks is skipped
    assert df.columns == ["#CHROM", "sample"]
    assert len(df) == len(datagen.WOMBAT_ROWS)


def test_read_cohort_wombat_columns_of_a_borrowed_header(data_dir: Path):
    headerless = datagen.build_headerless_cohort_wombat(data_dir)
    family = datagen.family_dir(data_dir, "1001") / "wombat"
    df = read_cohort_wombat(
        headerless,
        [family / "1001.rare.vep1.annotated.strict.tsv"],
        columns=["VEP_SYMBOL", "POS"],
        infer_schema_length=0,
    )
    assert df.columns == ["POS", "VEP_SYMBOL"]
    assert df["VEP_SYMBOL"].to_list() == ["SHANK3", "KIAA1109"]


def test_every_kind_is_indexed_from_the_cohort_tables(indexed: Path):
    entries = _manifest(indexed).values()
    kinds = sorted((e["kind"], e["cohort"], e["label"]) for e in entries)
    assert kinds == [
        ("roh", "DEMO", f"{datagen.BUNDLE}_{datagen.TAGS[0]}"),
        ("snv", "DEMO", "default"),
        ("snv", "DEMO", "strict"),  # the headerless one, its header borrowed
        ("snv", "OTHER", "default"),
        ("str", "DEMO", datagen.EH_LABEL),
        ("sv", "DEMO", "wisecondorx"),
    ]
    assert not any(e.get("error") for e in entries)


def test_snv_rows_one_per_variant_sample_gene_and_config(indexed: Path):
    df = gene_index.lookup_gene("snv", ["SHANK3"])
    rows = sorted(
        df.select("cohort", "config", "family_id", "sample", "variant").rows()
    )
    assert rows == [
        ("DEMO", "default", "1001", "S3", SHANK3_DE_NOVO),
        ("DEMO", "strict", "1001", "S3", SHANK3_DE_NOVO),
        ("OTHER", "default", "2001", "S9", SHANK3_DE_NOVO),
        ("OTHER", "default", "2001", "S9", "chr22:50675000:C:T"),
    ]
    assert df["data_dir"].unique().to_list() == [str(indexed)]


def test_sv_rows_carry_the_key_curation_uses(indexed: Path):
    [row] = gene_index.lookup_gene("sv", ["SCN2A"]).to_dicts()
    assert row["variant"] == "chr2:165114001-165316000:del"
    assert row["exonic"] is True
    assert row["family_id"] == "1001"


def test_str_and_roh_rows(indexed: Path):
    htt = gene_index.lookup_gene("str", ["HTT"])
    assert len(htt) == len(datagen.DEMO_PEDIGREE) - 1  # 1002 has no STR results
    assert htt.filter(pl.col("sample") == "S3")["variant"].item() == (
        "chr4:3074876-3074933:HTT"
    )
    roh = gene_index.lookup_gene("roh", ["genea"])
    assert sorted(roh["sample"].to_list()) == ["S3", "S6"]
    assert roh["chrom"].unique().to_list() == ["chr2"]


# ---------------------------------------------------------------------------
# Keeping it current
# ---------------------------------------------------------------------------


def test_refresh_reads_only_what_changed(indexed: Path):
    root = gene_index._store_root(str(indexed))
    files = {p: e["file"] for p, e in _manifest(indexed).items()}
    stamps = {p: (root / f).stat().st_mtime_ns for p, f in files.items()}

    other = indexed / "cohorts/OTHER/wombat/OTHER.rare.vep1.default.results.tsv"
    datagen.write_tsv(other, datagen.WOMBAT_COLUMNS, datagen.OTHER_WOMBAT_ROWS[:1])
    os.utime(other, (1, 1))
    wcx = indexed / "cohorts/DEMO/svs/wisecondorx/DEMO_aberrations.bed"
    wcx.unlink()
    gene_index.refresh_all(wait=True)

    after = _manifest(indexed)
    assert str(wcx) not in after
    assert not (root / files[str(wcx)]).exists()
    for path, stamp in stamps.items():
        if path in (str(other), str(wcx)):
            continue
        assert (root / after[path]["file"]).stat().st_mtime_ns == stamp
    assert after[str(other)]["rows"] == 1


def test_a_restart_serves_the_index_left_on_disk(indexed: Path):
    gene_index.clear_gene_index()
    assert gene_index.lookup_gene("snv", ["SHANK3"]).is_empty()
    gene_index.refresh_all(wait=True)  # nothing changed: nothing is read
    assert len(gene_index.lookup_gene("snv", ["SHANK3"])) == 4


def test_an_unreadable_table_is_reported_and_the_rest_kept(indexed: Path):
    bad = indexed / "cohorts/OTHER/wombat/OTHER.rare.vep1.broken.results.tsv"
    bad.write_text("\n1\t2\t3\n")  # headerless, and no family file lends one
    gene_index.refresh_all(wait=True)
    status = gene_index.index_status()
    assert "no header line" in status.errors[str(bad)]
    assert not status.building
    assert len(gene_index.lookup_gene("snv", ["SHANK3"])) == 4


def test_a_schema_change_rebuilds_everything(indexed: Path, monkeypatch):
    monkeypatch.setattr(gene_index, "INDEX_VERSION", gene_index.INDEX_VERSION + 1)
    gene_index.clear_gene_index()
    root = gene_index._store_root(str(indexed))
    assert gene_index._load_manifest(root) == {}


@pytest.mark.parametrize(
    ("text", "label"),
    [
        ("chr22:50674000:G:A", SHANK3_DE_NOVO),
        ("22-50674000-g-a", SHANK3_DE_NOVO),
        ("chrX_70000000_T_C", "chrX:70000000:T:C"),
        ("MT:73:A:G", "chrM:73:A:G"),
        ("SHANK3", None),
        ("chr22:50674000", None),
    ],
)
def test_parse_variant(text, label):
    parsed = gene_index.parse_variant(text)
    assert (parsed.label if parsed else None) == label


def test_lookup_variant(indexed: Path):
    df = gene_index.lookup_variant(gene_index.parse_variant("22-50674000-G-A"))
    assert sorted(set(df["sample"].to_list())) == ["S3", "S9"]


# ---------------------------------------------------------------------------
# Page rows
# ---------------------------------------------------------------------------


def test_gene_keys_follow_previous_symbols(snapshot):
    current, keys = gene_rows.gene_keys("kiaa1109", snapshot)
    assert current == "BLTP1"
    assert {"BLTP1", "KIAA1109"} <= keys


def test_snv_rows_join_configs_and_curation(indexed: Path):
    dirs = gene_rows.data_dir_infos()
    curation = gene_rows.load_curation(dirs)
    carriers = gene_rows.snv_carrier_rows(
        gene_index.lookup_gene("snv", ["SHANK3"]), curation, dirs, _scorer()
    )
    by_sample = {(r["Sample"], r["Variant"]): r for r in carriers}
    s3 = by_sample[("S3", SHANK3_DE_NOVO)]
    assert s3["Configs"] == "default, strict"
    assert s3["Validation"] == "present"
    assert s3["Diagnostic"] == "pathogenic"
    assert s3["_next"] == "/cohort/DEMO/family/1001?tab=wombat&wcfg=default"
    # The ignored "present" does not count
    assert by_sample[("S9", SHANK3_DE_NOVO)]["Validation"] == "absent"

    [shared, other] = gene_rows.snv_variant_rows(carriers)
    assert shared["Variant"] == SHANK3_DE_NOVO
    assert (shared["Carriers"], shared["Families"]) == (2, 2)
    assert shared["Cohorts"] == "DEMO, OTHER"
    assert shared["Validations"] == "absent 1 · present 1"
    assert other["Validations"] == "TODO 1"
    assert gene_rows.seen_summary(carriers) == {
        "individuals": 2,
        "families": 2,
        "cohorts": 2,
        "data_dirs": 1,
    }


def test_str_rows_keep_genotypes_past_normal(indexed: Path):
    dirs = gene_rows.data_dir_infos()
    frame = gene_index.lookup_gene("str", ["HTT"])
    rows = gene_rows.str_rows(frame, gene_rows.load_curation(dirs), dirs, _scorer())
    assert [(r["Sample"], r["Alleles"]) for r in rows] == [("S3", "18/40")]
    assert gene_rows.str_genotyped(frame) == 5


def test_diagnostics_found_by_gene_including_an_sv_overlap(indexed: Path):
    dirs = gene_rows.data_dir_infos()
    rows = gene_rows.diagnostic_rows(dirs, dirs, _scorer(), gene_keys={"SCN2A"})
    assert [(r["Type"], r["Diagnostic"], r["Cohort"]) for r in rows] == [
        ("SV", "uncertain", "DEMO")
    ]
    assert rows[0]["_next"] == "/family/1001"


def test_variant_validation_history_keeps_ignored_records(indexed: Path):
    dirs = gene_rows.data_dir_infos()
    rows = gene_rows.validation_history(dirs, dirs, "SNV", SHANK3_DE_NOVO)
    assert [(r["Sample"], r["Validation"], r["Ignored"]) for r in rows] == [
        ("S9", "absent", ""),
        ("S9", "present", "yes"),
        ("S3", "present", ""),
    ]


def test_rows_of_another_data_directory_link_through_it(app_config, data_dir, tmp_path):
    second = datagen.build_second_data_dir(tmp_path)
    data.init_all_data_stores(
        [DataDirectoryConfig(path=str(data_dir)), DataDirectoryConfig(path=str(second))]
    )
    gene_index.refresh_all(wait=True)
    dirs = gene_rows.data_dir_infos()
    carriers = gene_rows.snv_carrier_rows(
        gene_index.lookup_gene("snv", ["SHANK3"]),
        gene_rows.load_curation(dirs),
        dirs,
        _scorer(),
    )
    extra = next(r for r in carriers if r["Cohort"] == "EXTRA")
    assert extra["_dir"] == "1"
    assert extra["Data dir"] == "data2"
    assert gene_rows.seen_summary(carriers)["data_dirs"] == 2
