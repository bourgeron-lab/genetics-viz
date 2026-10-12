"""Small pure utilities: chromosomes, VEP, ClinVar, scores, loci, CNVs, labels, genes."""

import polars as pl
import pytest

from genetics_viz.utils import (
    clinvar,
    column_names,
    cytobands,
    gene_scoring,
    locus,
    pedigree_labels,
    score_colors,
    vep,
    view_presets,
    wisecondorx,
)


@pytest.mark.parametrize(
    "raw, norm",
    [
        ("chr1", "1"),
        ("1", "1"),
        ("Chr1", "1"),
        ("CHRX", "X"),
        ("chrM", "MT"),
        ("m", "MT"),
    ],
)
def test_norm_chrom(raw, norm):
    assert cytobands.norm_chrom(raw) == norm


def test_cytoband_range():
    assert cytobands.get_cytoband_range("chr1", 0, 1) == "1p36.33"
    spanning = cytobands.get_cytoband_range("Chr1", 0, 7_000_000)
    assert spanning.startswith("1p36.33-") and spanning != "1p36.33"
    assert cytobands.get_cytoband_range("chrZ", 0, 1) == ""


def test_vep_priorities():
    assert vep.get_consequence_impact("stop_gained") == "HIGH"
    worst = vep.get_highest_consequence_term(
        "missense_variant&splice_region_variant,synonymous_variant"
    )
    assert worst == "missense_variant"
    assert vep.get_highest_priority_consequence("") == 9999
    assert vep.format_consequence_display("splice_region_variant") == "splice region"


def test_clinvar_colors():
    known = next(iter(clinvar.CLINVAR_COLORS))
    assert clinvar.get_clinvar_color(known.upper()) == clinvar.get_clinvar_color(known)
    assert clinvar.get_clinvar_color("") == "#757575"


def test_score_colors():
    low = score_colors.get_score_color("REVEL", -1)
    high = score_colors.get_score_color("revel", 10)
    assert low["color"].lower() == "#4caf50" and low != high
    assert score_colors.get_score_color("not_a_score", 1) is None
    assert score_colors.interpolate_color("#000000", "#ffffff", 0.5).lower() in (
        "#7f7f7f",
        "#808080",
    )


def test_genomic_sort_key_orders_chromosomes_naturally():
    keys = ["chr10:5:A:T", "chrX:1", "chr2:9", "chrM:3", "1:7", None, "garbage"]
    ordered = sorted(keys, key=column_names.genomic_sort_key)
    assert ordered[:5] == ["1:7", "chr2:9", "chr10:5:A:T", "chrX:1", "chrM:3"]


def test_view_preset_selection():
    presets = [
        {"name": "A", "keywords": ["denovo"]},
        {"name": "B", "keywords": ["homalt", "HOM"]},
    ]
    assert view_presets.select_preset_for_config("rare_HOMALT", presets)["name"] == "B"
    assert view_presets.select_preset_for_config("other", presets)["name"] == "A"
    assert view_presets.select_preset_for_config("x", [])["name"] == "Default"


@pytest.mark.parametrize(
    "query, kind",
    [
        ("chr1:100:a:gc", "exact_variant"),
        ("1:100-200", "range"),
        ("chrX:100", "exact_position"),
        ("ensg00000164099", "gene_id"),
        ("SHANK*", "gene_wildcard"),
        ("shank3", "gene_name"),
    ],
)
def test_parse_locus_query(query, kind):
    parsed = locus.parse_locus_query(query)
    assert parsed["type"] == kind
    if kind == "exact_variant":
        assert (parsed["chrom"], parsed["ref"], parsed["alt"]) == ("chr1", "A", "GC")
    if kind == "range":
        assert parsed["chrom"] == "chr1" and (parsed["start"], parsed["end"]) == (
            100,
            200,
        )


def test_cnv_classification_and_sv_type():
    cfg = wisecondorx.WISECONDORX_CONFIG
    rl, rg = cfg["robust_loss"], cfg["robust_gain"]
    assert (
        wisecondorx.classify_cnv(rl["ratio_threshold"] - 1, rl["zscore_threshold"] - 1)
        == rl["label"]
    )
    assert (
        wisecondorx.classify_cnv(rg["ratio_threshold"] + 1, rg["zscore_threshold"] + 1)
        == rg["label"]
    )
    assert wisecondorx.classify_cnv(0, 0) == "Below threshold"
    assert wisecondorx.classify_cnv("x", 1) == "N/A"
    assert wisecondorx.infer_sv_type({"call": "Robust GAIN"}) == "dup"
    assert wisecondorx.infer_sv_type({"type": "del"}) == "del"
    assert wisecondorx.infer_sv_type({"ratio": "0.4"}) == "dup"
    assert wisecondorx.infer_sv_type({}) == "del"


def test_wisecondorx_bed(tmp_path):
    path = tmp_path / "x.bed"
    path.write_text(
        "chr\tstart\tend\tratio\tzscore\ttype\tbarcode\n1\t100\t200\t-1.2\t-9\tloss\tS1\nshort\trow\n"
    )
    df = wisecondorx.parse_wisecondorx_bed(path)
    assert df is not None and len(df) == 1
    row = df.row(0, named=True)
    assert row["sample"] == "S1" and str(row["chr"]).startswith("chr")


@pytest.mark.parametrize(
    "raw, norm",
    [("1.0", "1"), ("2", "2"), ("-9.0", "-9"), ("0", "-9"), (None, "-9"), ("-", "-9")],
)
def test_pedigree_codes(raw, norm):
    assert pedigree_labels.normalize_pheno(raw) == norm
    assert pedigree_labels.normalize_sex(raw) == norm


def test_every_pedigree_code_has_a_label_and_colour():
    for code in pedigree_labels.PHENO_ORDER:
        assert (
            code in pedigree_labels.PHENO_LABELS
            and code in pedigree_labels.PHENO_COLORS
        )
    for code in pedigree_labels.SEX_ORDER:
        assert code in pedigree_labels.SEX_LABELS and code in pedigree_labels.SEX_COLORS


def test_gene_scoring_colours_and_lookup():
    scorer = gene_scoring.get_gene_scorer()
    assert gene_scoring.GeneScorer.score_to_color(0) == "#ffffff"
    assert gene_scoring.GeneScorer.score_to_color(
        100
    ) == gene_scoring.GeneScorer.score_to_color(8)
    assert scorer.get_gene_score_and_sets("NOT_A_GENE_XYZ") == (0.0, [])
    assert gene_scoring.get_gene_scorer() is scorer


def test_schema_overrides_use_polars_types():
    for dtype in column_names.get_schema_overrides().values():
        assert isinstance(dtype, pl.DataType) or (
            isinstance(dtype, type) and issubclass(dtype, pl.DataType)
        )
