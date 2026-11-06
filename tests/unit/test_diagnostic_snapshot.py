"""Variant details saved with a diagnostic (utils/diagnostic_snapshot.py), the
transcript picked before Wombat rows are merged, and the file upgrade."""

import csv

import datagen
import polars as pl
import pytest

from genetics_viz.components import diagnostic_loader as dl
from genetics_viz.utils.diagnostic_snapshot import (
    SNAPSHOT_COLUMNS,
    has_snapshot,
    snv_snapshot,
    vedette_snapshot,
    wisecondorx_snapshot,
)
from genetics_viz.utils.wombat import aggregate_by_variant

PHENOTYPES = {"S1": "1", "S2": "1", "S3": "2"}


@pytest.fixture(autouse=True)
def curator(monkeypatch):
    monkeypatch.setattr(dl, "can_write", lambda: True)


@pytest.fixture
def family_rows(data_dir):
    path = (
        datagen.family_dir(data_dir, "1001")
        / "wombat"
        / "1001.rare.vep1.annotated.default.tsv"
    )
    df = pl.read_csv(path, separator="\t", infer_schema=False)
    return aggregate_by_variant(df).to_dicts()


def _row(rows, pos):
    return next(r for r in rows if r["POS"] == str(pos))


def test_the_mane_transcript_is_picked_before_rows_merge(family_rows):
    shank3 = _row(family_rows, 50674000)
    assert shank3["n_grouped"] == 2
    # Merged, the two transcripts' HGVS are a list nobody can pair up...
    assert shank3["VEP_HGVSc"].count(",") == 1
    # ...so the MANE Select one is taken from its own row first.
    assert shank3["_hgvs_transcript"] == "ENST00000710353.1"
    assert shank3["_hgvs_mane"] == "NM_033517.1"
    assert shank3["_hgvsc"] == "ENST00000710353.1:c.3679G>A"
    assert shank3["_hgvsp"] == "ENSP00000518234.1:p.Gly1227Arg"


def test_canonical_is_the_fallback_and_nothing_else_is_picked():
    df = pl.DataFrame(
        {
            "#CHROM": ["1", "1", "2"], "POS": [1, 1, 5], "REF": ["A"] * 3,
            "ALT": ["C"] * 3, "sample": ["S"] * 3,
            "VEP_Feature": ["ENST1", "ENST2", "ENST9"],
            "VEP_HGVSc": ["ENST1:c.1A>C", "ENST2:c.5A>C", "ENST9:c.9A>C"],
            "VEP_MANE_SELECT": [".", ".", "."],
            "VEP_CANONICAL": [".", "YES", "."],
        }
    )  # fmt: skip
    first, second = aggregate_by_variant(df).to_dicts()
    assert first["_hgvsc"] == "ENST2:c.5A>C"
    assert first["_hgvs_mane"] is None
    assert second["_hgvsc"] is None
    assert "_hgvsc" not in aggregate_by_variant(df.drop("VEP_HGVSc")).columns


def test_snv_snapshot_of_a_de_novo(family_rows):
    shank3 = {
        **_row(family_rows, 50674000),
        "MOI": "AD ✓",
        "MoiBadges": [{"label": "AD ✓", "tooltip": "SHANK3: consistent\nmore"}],
    }
    snap = snv_snapshot(
        shank3, source="Wombat default", rows=family_rows, phenotypes=PHENOTYPES
    )
    assert set(snap) == set(SNAPSHOT_COLUMNS)
    assert snap["Source"] == "Wombat default"
    assert (snap["MANE"], snap["HGVSc"], snap["HGVSp"]) == (
        "NM_033517.1",
        "c.3679G>A",
        "p.Gly1227Arg",
    )
    assert (snap["Zygosity"], snap["Segregation"], snap["Carriers"]) == (
        "het",
        "de novo",
        "",
    )
    assert snap["ClinVar"] == "pathogenic"
    # No AF, but a FAF95: in gnomAD, so not "absent".
    assert (snap["gnomAD_AF"], snap["gnomAD_FAF95"]) == ("", "0.0")
    assert snap["MOI"] == "SHANK3: consistent"
    no_gnomad = {**shank3, "fafmax_faf95_max_genomes": None}
    assert snv_snapshot(no_gnomad, source="w")["gnomAD_AF"] == "absent"
    # Without badges, the cell as it is.
    assert snv_snapshot({**shank3, "MoiBadges": []}, source="w")["MOI"] == "AD ✓"


def test_snv_snapshot_of_inherited_and_hemizygous_calls(family_rows):
    paternal = snv_snapshot(
        _row(family_rows, 32000000), source="w", rows=family_rows,
        phenotypes=PHENOTYPES,
    )  # fmt: skip
    assert paternal["Segregation"] == "paternal"
    assert paternal["Carriers"] == "S1 (unaffected, 0/1)"
    assert paternal["HGVSc"] == ""  # no MANE, no canonical: none picked

    plain = snv_snapshot(_row(family_rows, 1000), source="w")
    assert plain["HGVSp"] == "p.Leu4="  # VEP writes "=" as %3D
    assert plain["gnomAD_AF"] == "1.2e-05"

    xlinked = snv_snapshot(_row(family_rows, 70000000), source="w")
    assert (xlinked["Zygosity"], xlinked["Segregation"]) == ("hemi", "maternal")


def test_a_row_with_one_transcript_needs_no_pick():
    row = {
        "VEP_HGVSc": "ENST5.1:c.2del",
        "VEP_HGVSp": "ENSP5.1:p.Ter1?",
        "VEP_MANE_SELECT": "NM_5.2",
    }
    snap = snv_snapshot(row, source="file")
    assert (snap["Transcript"], snap["MANE"], snap["HGVSc"]) == (
        "ENST5.1",
        "NM_5.2",
        "c.2del",
    )
    two = snv_snapshot({"VEP_HGVSc": "ENST5.1:c.2del,ENST6.1:c.9del"}, source="f")
    assert two["HGVSc"] == ""


def test_wisecondorx_snapshot():
    row = {
        "gene": "SCN2A:exonic,SCN3A:intronic",
        "ratio": -1.0,
        "zscore": -8.0,
        "_original_locus": "chr2:165114001-165316000",
    }
    snap = wisecondorx_snapshot(row)
    assert snap["Source"] == "WisecondorX"
    assert snap["GeneOverlap"] == "SCN2A:exonic, SCN3A:intronic"
    assert snap["CopyNumber"] == "≈1.0 (log2 ratio -1.00, z -8.0)"
    # On X the expected ploidy depends on the karyotype: the ratio only.
    on_x = wisecondorx_snapshot({**row, "_original_locus": "chrX:1-2"})
    assert on_x["CopyNumber"] == "log2 ratio -1.00, z -8.0"


def test_vedette_snapshot():
    base = {"site_id": "7", "genes": "RERE", "inheritance": "de novo"}
    rows = [
        {**base, "sample": "S3", "GT": "0/1", "CN": 1},
        {**base, "sample": "S1", "GT": "0/1", "CN": 1},
        {**base, "site_id": "8", "sample": "S2", "GT": "1/1", "CN": 0},
    ]
    snap = vedette_snapshot(rows[0], source="Vedette 0.4", rows=rows)
    assert (snap["Zygosity"], snap["Segregation"], snap["CopyNumber"]) == (
        "het",
        "de novo",
        "1",
    )
    assert snap["Carriers"] == "S1 (phenotype unknown, 0/1)"
    assert vedette_snapshot(rows[2], source="v")["Zygosity"] == "hom"
    assert vedette_snapshot({**rows[2], "CN": 0}, source="v")["CopyNumber"] == "0"


def test_saving_upgrades_an_older_file(tmp_path):
    path = tmp_path / "diagnostics" / "snvs.tsv"
    datagen.write_tsv(
        path,
        datagen.DIAGNOSTIC_COLUMNS,
        [["1001", "chr1:1:A:T", "G", "x", "S3", "u", "t0", "", "0", "benign"]],
    )
    details = {"Source": "Wombat default", "HGVSc": "c.1A>T", "Zygosity": "het"}
    dl.save_diagnostic_entry(
        path, "1001", "chr1:1:A:T", "G", "x", "S3", "u", "t1", "", "pathogenic",
        details=details,
    )  # fmt: skip
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        old, new = list(reader)
    assert reader.fieldnames == dl.DIAGNOSTIC_HEADER
    assert not has_snapshot(old) and old["Diagnostic"] == "benign"
    assert (new["Diagnostic"], new["HGVSc"], new["Zygosity"]) == (
        "pathogenic",
        "c.1A>T",
        "het",
    )
    assert None not in new  # no field beyond the header

    [entry] = dl.load_family_diagnostics(tmp_path, "1001", ["S3"])
    assert entry["Diagnostic"] == "conflicting"
    assert [has_snapshot(r) for r in entry["_entries"]] == [True, False]


def test_mane_is_not_given_to_another_transcripts_hgvs():
    """#37: the MANE transcript has no c.; the one c. left is another's."""
    row = {
        "VEP_Feature": "ENST1.1,ENST2.1",
        "VEP_HGVSc": "ENST2.1:c.-5A>G",
        "VEP_MANE_SELECT": "NM_1.1",
    }
    snap = snv_snapshot(row, source="w")
    assert (snap["Transcript"], snap["HGVSc"], snap["MANE"]) == (
        "ENST2.1",
        "c.-5A>G",
        "",
    )
    one = snv_snapshot({**row, "VEP_Feature": "ENST2.1"}, source="w")
    assert one["MANE"] == "NM_1.1"


def test_vedette_snapshot_words_codes_and_tags_exonic_genes():
    row = {
        "sample": "S3", "GT": "0/1", "CN": 1, "site_id": "7",
        "inheritance": "de novo (low conf.)", "_inheritance": "denovo_lowconf",
        "genes": "RERE, ENO1", "_genes": ["RERE", "ENO1"], "_exonic": ["RERE"],
    }  # fmt: skip
    snap = vedette_snapshot(row, source="v")
    assert snap["Segregation"] == "de novo (low confidence)"
    assert snap["GeneOverlap"] == "RERE:exonic, ENO1:genic"
    # A founder's "none" is no segregation, not a word in the report
    founder = vedette_snapshot({**row, "_inheritance": "none"}, source="v")
    assert founder["Segregation"] == ""
    assert (
        vedette_snapshot({**row, "_inheritance": "biparental"}, source="v")[
            "Segregation"
        ]
        == "both"
    )


def test_wisecondorx_gene_overlap_is_stable():
    from genetics_viz.utils.wisecondorx import gene_overlap

    assert gene_overlap("B,A,C", "C,A") == "A:exonic,C:exonic,B:genic"
    assert gene_overlap(None, "") == ""
