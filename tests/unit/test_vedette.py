"""Vedette calls: the family export (v1.4 format and before), carrier rows,
filters and presets, link, own call.

Also the read evidence Vedette's small calls need: deletions inside reads.
"""

import gzip

import pysam

from datagen import (
    CRAM_INREAD_DELS,
    VEDETTE_CALLER,
    VEDETTE_DEL,
    VEDETTE_NONPASS,
    VEDETTE_SITE_KEY,
    VEDETTE_START,
    cram_reference,
    family_dir,
    sample_dir,
)
from genetics_viz.components.validation_loader import (
    VEDETTE_VALIDATION_HEADER,
    load_validation_map,
    merge_validation_maps,
)
from genetics_viz.utils import vedette as v
from genetics_viz.utils.sv_breakpoints import suggest_breakpoints
from genetics_viz.utils.sv_read_plots import SERIES_NAMES, read_option
from genetics_viz.utils.sv_reads import load_read_evidence
from genetics_viz.utils.sv_view import read_view_windows
from genetics_viz.utils.tsv import append_tsv_rows
from genetics_viz.utils.validation_files import VALIDATION_TYPES, validation_path

MOSAIC = "mosaic_S3_chr2_180000001"
SUB = "sub_S3_chr2_175000001"


def _rows(data_dir):
    """As the tab reads them: the Parquet file when it is there."""
    info = v.read_export_info(data_dir, "1001")
    table, members, details = v.read_family_rows(data_dir, "1001", info)
    return v.carrier_rows(
        table,
        labels={"S1": "father", "S2": "mother"},
        members=members,
        own_details=details,
    )


def _by_id(rows):
    return {(r["site_id"], r["sample"]): r for r in rows}


# ---------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------


def test_impact_label_of_older_exports_is_split():
    assert v.parse_impact("+0.90 — élevé · DEL exonique RERE (NDD, s_het 0.1)") == (
        0.9,
        "high",
        "DEL exonique RERE (NDD, s_het 0.1)",
    )
    assert v.parse_impact("-0.60 — aucun gène codant") == (-0.6, "none", "")
    assert v.parse_impact("+0.00 — faible · intronique")[1] == "low"
    assert v.parse_impact("-1.00 — bénin · région ClinGen")[1] == "benign"
    # A level Vedette adds later is kept as written, coloured grey.
    assert v.parse_impact("+0.20 — nouveau · x")[1] == "nouveau"
    assert v.level_info("nouveau")["color"] == v.UNKNOWN_COLOR
    assert v.parse_impact(None) == (None, "", "")


def test_codes_are_worded():
    assert v.parse_inheritance("B1:biparental; B2:denovo_lowconf") == {
        "B1": "biparental",
        "B2": "denovo_lowconf",
    }
    assert v.parse_inheritance("") == {}
    assert v.split_site_genes("SKI, MORN1, PANK4 (+5)") == (
        ["SKI", "MORN1", "PANK4"],
        5,
    )
    assert v.split_site_genes(None) == ([], 0)
    assert v.inheritance_info("karyotype_XXY")["label"] == "karyotype XXY"
    assert v.inheritance_info("other")["color"] == v.UNKNOWN_COLOR
    assert v.kind_info("subthreshold")["label"] == "sub-threshold"
    assert v.own_call_info("none")["label"] == "re-genotyped only"
    assert v.level_info("benign")["label"] == "Benign (ClinGen)"


def test_one_row_per_event_and_carrier(data_dir):
    rows = _rows(data_dir)
    assert [(r["site_id"], r["sample"]) for r in rows] == [
        ("vedette_site_1", "S3"),
        ("vedette_site_2", "S1"),
        ("vedette_site_2", "S3"),
        ("vedette_site_4", "S1"),
        ("vedette_site_4", "S2"),
        ("vedette_site_4", "S3"),
        ("vedette_site_5", "S2"),
        ("vedette_site_5", "S3"),
        (SUB, "S3"),
        (MOSAIC, "S3"),
        (VEDETTE_NONPASS, "S1"),
        (VEDETTE_NONPASS, "S3"),
    ]
    by = _by_id(rows)
    first = by[("vedette_site_1", "S3")]
    # 1-based inclusive bounds, as the export writes them.
    assert (first["start"], first["end"], first["size"]) == (
        VEDETTE_START,
        VEDETTE_DEL[1],
        120,
    )
    # Keyed by the cohort site (site_key), not by the row's refined bounds.
    assert first["_variant_key"] == VEDETTE_SITE_KEY != first["_locus_key"]
    assert first["_locus_key"] == f"chr2:{VEDETTE_START}-{VEDETTE_DEL[1]}:del"
    # The cohort site's INFO, and the carrier's own call (v1.9.0).
    assert (first["af"], first["nsamples"], first["precise_frac"]) == (0.002, 3, 1.0)
    assert (first["SQ"], first["source"], first["_precise"]) == (12.5, "JUNCTION", True)
    assert first["_own_call_is_member"] and first["_own_detail"]["pe"] == 3
    assert (first["GT"], first["CN"], first["GQ"], first["DNQ"]) == ("0/1", 1, 45, 38)
    assert first["inheritance"] == "de novo" and first["impact_level"] == "High"
    assert (first["impact_score"], first["impact_reason"]) == (
        0.9,
        "DEL exonique SCN2A (NDD candidat, s_het 0.12)",
    )
    assert (first["kind"], first["own_call"], first["support"]) == (
        "site",
        "own call PASS",
        "CN agrees",
    )
    assert first["_exonic"] == ["SCN2A"] and first["cohort_carrier_frac"] == 0.004
    # Every member's genotype; the father named, the child's DNQ in the tooltip.
    assert [b["label"] for b in first["family_badges"]] == ["0/0·2", "0/0·2", "0/1·1"]
    assert first["family_badges"][0]["tooltip"].startswith("S1 (father): GT 0/0")
    assert "DNQ 38" in first["family_badges"][2]["tooltip"]
    # Two cohort sites collapsed into one event; a founder has no inheritance.
    dup = by[("vedette_site_2", "S1")]
    assert dup["site_ids"] == ["vedette_site_2", "vedette_site_3"]
    assert dup["n_sites"] == 2 and dup["inheritance"] == ""
    assert dup["_genes"] == ["SCN2A"] and dup["_exonic"] == []
    # The complete gene list, where the old export wrote "(+1)".
    xrow = by[("vedette_site_5", "S3")]
    assert len(xrow["_genes"]) == 7 and "(+" not in xrow["genes"]
    sub = by[(SUB, "S3")]
    assert (sub["kind"], sub["filter"], sub["GQ"]) == ("sub-threshold", "LowSQ", 15)
    assert sub["cohort_carrier_frac"] is None
    # A mosaic row: carried by S3 though no genotype is written.
    mosaic = by[(MOSAIC, "S3")]
    assert (mosaic["GT"], mosaic["GQ"], mosaic["mosaic_cn"]) == ("", None, 2.4)
    assert (
        mosaic["impact_score"] is None and mosaic["family_badges"][0]["label"] == "-·-"
    )


def test_export_info_and_member_qc(data_dir):
    info = v.read_export_info(data_dir, "1001")
    assert (info.schema_version, info.caller, info.model_version) == (
        1,
        VEDETTE_CALLER,
        "v5n",
    )
    assert info.members == ["S1", "S2", "S3"] and not info.newer_schema
    assert v.caller_version(data_dir, "1001") == VEDETTE_CALLER
    qc = {m.sample: m for m in v.read_member_qc(data_dir, "1001")}
    assert qc["S3"].flags == ["carrier_load", "denovo_outlier", "wavy"]
    assert (qc["S3"].carrier_excess, qc["S3"].denovo_candidates) == (5.47, 458)
    assert qc["S1"].flags == [] and qc["S1"].denovo_candidates is None
    assert v.qc_flag_label("wavy") == "wavy coverage"
    # A family without an export.
    assert v.read_export_info(data_dir, "1002").caller == ""
    assert v.read_member_qc(data_dir, "1002") == []


def test_parquet_and_tsv_give_the_same_rows(data_dir):
    info = v.read_export_info(data_dir, "1001")
    assert (info.pon, info.parquet, info.nonpass_max_carrier_frac) == (
        "dp_founders300",
        "1001.vedette.parquet",
        0.2,
    )
    table, members, details = v.read_family_rows(data_dir, "1001", info)
    assert isinstance(table, list) and members == ["S1", "S2", "S3"]
    assert set(details) == {("vedette_site_1", "S3"), ("vedette_site_5", "S3")}
    from_parquet = v.carrier_rows(table, members=members, own_details=details)
    from_tsv = v.carrier_rows(
        v.read_vedette_table(v.family_vedette_file(data_dir, "1001")),
        members=info.members,
    )

    def shown(rows):
        return [{k: val for k, val in r.items() if k != "_own_detail"} for r in rows]

    # Same rows; only the Parquet has the own calls' details.
    assert shown(from_parquet) == shown(from_tsv)
    # Without the Parquet file, the TSV is read.
    v.family_parquet_file(data_dir, "1001", info).unlink()
    table, _, details = v.read_family_rows(data_dir, "1001", info)
    assert not isinstance(table, list) and details == {}


def test_own_call_comes_from_the_export(data_dir):
    by = _by_id(_rows(data_dir))
    # Parquet: the whole own call, no VCF opened.
    own = v.resolve_own_call(data_dir / "nowhere", by[("vedette_site_1", "S3")])
    assert (own.origin, own.call_id, own.start, own.end) == (
        "parquet",
        "vedette_DEL_7",
        VEDETTE_START,
        VEDETTE_DEL[1],
    )
    assert own.passed and own.fields["CG"] == 4 and own.fields["SQ"] == 12.5
    filtered = v.resolve_own_call(data_dir / "nowhere", by[("vedette_site_5", "S3")])
    assert filtered.filters == ["LowSQ"] and filtered.cipos == (-250, 250)
    assert not filtered.precise
    # "none": re-genotyped only, nothing to look up.
    assert v.resolve_own_call(data_dir, by[("vedette_site_2", "S3")]) is None
    # Precision from the export, else the own call, else the size.
    assert v.row_is_precise(by[("vedette_site_5", "S3")], None) is False
    assert v.row_is_precise(by[(MOSAIC, "S3")], None) is False  # 2 Mb
    # From the TSV alone: the code says pass, the VCF gives the rest.
    info = v.read_export_info(data_dir, "1001")
    tsv = _by_id(
        v.carrier_rows(
            v.read_vedette_table(v.family_vedette_file(data_dir, "1001")),
            members=info.members,
        )
    )
    own = v.resolve_own_call(data_dir, tsv[("vedette_site_1", "S3")])
    assert (own.origin, own.call_id) == ("vcf", "vedette_DEL_7")
    # No VCF either: what the TSV says.
    partial = v.resolve_own_call(data_dir / "nowhere", tsv[("vedette_site_1", "S3")])
    assert (partial.origin, partial.call_id, partial.source) == ("tsv", "", "JUNCTION")
    assert partial.passed and partial.fields == {"SQ": 12.5}


def test_a_newer_schema_is_flagged():
    assert v.ExportInfo(schema_version=2).newer_schema
    assert not v.ExportInfo(schema_version=None).newer_schema


def test_exports_before_v1_4_still_open(tmp_path):
    """One row per PASS site, French impact label, truncated genes, an
    inheritance column, no JSON: the caller version from the family VCF."""
    data = tmp_path / "data"
    vdir = family_dir(data, "F1") / "svs" / "vedette"
    vdir.mkdir(parents=True)
    (vdir / "F1.vedette.tsv").write_text(
        "site_id\tcontig\tstart\tend\tsvtype\tsize\tfilter\tcohort_af\t"
        "cohort_carrier_frac\tcohort_ncarrier\tgenes\timpact\t"
        "GT_A\tCN_A\tGQ_A\tGT_B\tCN_B\tGQ_B\tinheritance\n"
        "vedette_site_9\tchr1\t1000\t2000\tDEL\t1001\tPASS\t0.002\t0.004\t3\t"
        "SKI, MORN1 (+4)\t+0.30 — faible · DEL exonique SKI\t"
        "0/1\t1\t30\t0/0\t2\t40\tA:denovo\n"
    )
    with gzip.open(vdir / "F1.vedette.vcf.gz", "wt") as f:
        f.write("##fileformat=VCFv4.2\n##source=vedette-cohort-1.1.0\n#CHROM\n")
    info = v.read_export_info(data, "F1")
    assert info.caller == "vedette-cohort-1.1.0" and info.schema_version is None
    [row] = v.carrier_rows(v.read_vedette_table(v.family_vedette_file(data, "F1")))
    assert (row["site_id"], row["sample"], row["inheritance"]) == (
        "vedette_site_9",
        "A",
        "de novo",
    )
    assert (row["impact_score"], row["_impact_level"]) == (0.3, "low")
    assert row["_genes"] == ["SKI", "MORN1"] and row["genes"].endswith("(+4)")
    assert (row["cohort_af"], row["cohort_ncarrier"], row["kind"]) == (0.002, 3, "site")


def test_missing_and_unknown_columns_degrade(tmp_path):
    path = tmp_path / "F.vedette.tsv"
    path.write_text(
        "row_id\tcontig\tstart\tend\tsvtype\tnew_col\tGT_A\tGT_B\n"
        "s1\tchr1\t100\t300\tDEL\tx\t0/1\t0/0\n"
    )
    [row] = v.carrier_rows(v.read_vedette_table(path))
    assert row["sample"] == "A" and row["size"] == 201
    assert row["new_col"] == "x" and row["impact_score"] is None
    assert row["CN"] is None and row["cohort_carrier_frac"] is None


# ---------------------------------------------------------------------------
# Filters, presets, link
# ---------------------------------------------------------------------------


def test_presets(data_dir):
    rows = _rows(data_dir)

    def keys(f, **kw):
        return [(r["site_id"], r["sample"]) for r in v.apply_filters(rows, f, **kw)]

    assert v.default_preset() == "rare"
    # Rare: carrier fraction <= 1 %, GQ >= 20. The sub-threshold call's GQ is
    # 15; the mosaic row has neither value, and a missing value filters nothing.
    assert keys(v.VedetteFilters.preset("rare")) == [
        ("vedette_site_1", "S3"),
        ("vedette_site_5", "S2"),
        ("vedette_site_5", "S3"),
        (MOSAIC, "S3"),
    ]
    assert keys(v.VedetteFilters.preset("denovo")) == [("vedette_site_1", "S3")]
    assert len(keys(v.VedetteFilters.preset("impact"))) == 3
    assert len(keys(v.VedetteFilters.preset("all"))) == 10
    # Non-PASS cohort sites are hidden unless asked for, whatever the preset.
    shown = v.VedetteFilters.preset("all")
    shown.show_nonpass = True
    assert len(keys(shown)) == 12
    impact = v.VedetteFilters.preset("impact")
    impact.show_nonpass = True
    assert (VEDETTE_NONPASS, "S3") in keys(impact)
    assert v.is_nonpass(_by_id(rows)[(VEDETTE_NONPASS, "S1")])
    assert not v.is_nonpass(_by_id(rows)[(SUB, "S3")])  # LowSQ is its kind
    # The page's member selection applies on top.
    assert keys(v.VedetteFilters.preset("all"), members=["S1"]) == [
        ("vedette_site_2", "S1"),
        ("vedette_site_4", "S1"),
    ]


def test_each_filter(data_dir):
    rows = _rows(data_dir)
    f = v.VedetteFilters.preset("all")
    f.types = ["DUP"]
    assert {r["site_id"] for r in v.apply_filters(rows, f)} == {
        "vedette_site_2",
        MOSAIC,
    }
    f = v.VedetteFilters.preset("all")
    f.kinds = ["subthreshold"]
    assert {r["site_id"] for r in v.apply_filters(rows, f)} == {SUB}
    f = v.VedetteFilters.preset("all")
    f.own_calls = ["none"]  # the carrier itself re-genotyped only (v1.9.0)
    assert {(r["site_id"], r["sample"]) for r in v.apply_filters(rows, f)} == {
        ("vedette_site_2", "S1"),
        ("vedette_site_2", "S3"),
        ("vedette_site_4", "S2"),
        ("vedette_site_4", "S3"),
        ("vedette_site_5", "S2"),
    }
    f.own_calls = ["filtered"]
    assert {(r["site_id"], r["sample"]) for r in v.apply_filters(rows, f)} == {
        ("vedette_site_4", "S1"),
        ("vedette_site_5", "S3"),
    }
    f = v.VedetteFilters.preset("all")
    f.min_size, f.max_size = 100, 350
    assert {r["site_id"] for r in v.apply_filters(rows, f)} == {
        "vedette_site_1",
        "vedette_site_5",
    }
    f = v.VedetteFilters.preset("all")
    f.inheritance = [""]  # "not given": founders, and rows without a trio code
    assert len(v.apply_filters(rows, f)) == 6
    f = v.VedetteFilters.preset("all")
    f.impact_levels = ["none"]
    assert {r["site_id"] for r in v.apply_filters(rows, f)} == {"vedette_site_4"}
    f = v.VedetteFilters.preset("all")
    assert len(v.apply_filters(rows, f, gene_symbols=frozenset({"SHANK3"}))) == 2


def test_link_round_trip():
    def through(f):
        params = v.filters_to_params(f)

        def get(key):
            value = params.get(key)
            if isinstance(value, list):
                return ",".join(map(str, value))
            return None if value is None else str(value)

        def get_list(key):
            value = params.get(key)
            return None if value is None else [str(x) for x in value]

        return params, v.filters_from_params(get, get_list)

    params, back = through(v.VedetteFilters.preset("denovo"))
    assert params["vpre"] == "denovo" and back.matches_preset("denovo")

    custom = v.VedetteFilters.preset("rare")
    custom.max_freq, custom.min_size, custom.inheritance = 0.05, 200, ["maternal", ""]
    custom.kinds, custom.own_calls = ["site", "mosaic"], ["pass"]
    custom.show_nonpass = True
    custom.genesets, custom.validations = ["tier:1"], ["present"]
    params, back = through(custom)
    assert params["vpre"] == "custom"
    assert (back.max_freq, back.min_size, back.max_size) == (0.05, 200, None)
    assert back.kinds == ["site", "mosaic"] and back.own_calls == ["pass"]
    assert back.show_nonpass and params["vnp"] == "1"
    _, plain = through(v.VedetteFilters.preset("rare"))
    assert not plain.show_nonpass
    assert back.inheritance == ["maternal", ""] and back.min_gq == 20
    assert back.genesets == ["tier:1"] and back.validations == ["present"]

    # No link: the opening preset; an unknown preset name falls back to it.
    assert v.filters_from_params(lambda k: None, lambda k: None).matches_preset("rare")
    assert v.filters_from_params(
        lambda k: "nope" if k == "vpre" else None, lambda k: None
    ).matches_preset("rare")


# ---------------------------------------------------------------------------
# Own call, validations
# ---------------------------------------------------------------------------


def test_own_call_is_matched_by_type_and_overlap(data_dir):
    vcf = v.sample_vedette_vcf(data_dir, "S3")
    own = v.find_own_call(vcf, "chr2", VEDETTE_START, VEDETTE_DEL[1], "DEL")
    assert own.call_id == "vedette_DEL_7" and own.passed and own.precise
    assert own.model == "v5n" and own.fields["CG"] == 4
    # ##vedette_coords=padded: POS is the base before, read back as the export's.
    assert (own.start, own.end) == (VEDETTE_START, VEDETTE_DEL[1])
    depth = v.find_own_call(vcf, "chrX", 2_800_001, 2_800_300, "DEL")
    assert (depth.start, depth.end) == (2_800_011, 2_800_290)
    assert not depth.precise and depth.filters == ["LowSQ"]
    assert depth.cipos == (-250, 250)
    # Another type, or no overlap: re-genotyped only.
    assert v.find_own_call(vcf, "chr2", *VEDETTE_DEL, "DUP") is None
    assert v.find_own_call(vcf, "chr2", 165_200_000, 165_200_400, "DUP") is None
    assert (
        v.find_own_call(v.sample_vedette_vcf(data_dir, "S1"), "chr2", 1, 2, "DEL")
        is None
    )
    assert v.assume_precise(300) and not v.assume_precise(5000)


def test_own_call_of_an_older_vcf(tmp_path):
    """Before v1.4.0: POS on the first affected base, END one past the last."""
    plain = tmp_path / "S.filtered.vcf"
    plain.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="e">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "chr1\t1001\told\tN\t<DEL>\t.\tPASS\tEND=1101;SVTYPE=DEL\tGT\t0/1\n"
    )
    pysam.tabix_compress(str(plain), f"{plain}.gz", force=True)
    pysam.tabix_index(f"{plain}.gz", preset="vcf", force=True)
    own = v.find_own_call(tmp_path / "S.filtered.vcf.gz", "chr1", 1001, 1100, "DEL")
    assert (own.start, own.end, own.overlap) == (1001, 1100, 1.0)


def test_vedette_validations_are_a_kind_of_their_own(data_dir):
    assert "Vedette" in VALIDATION_TYPES
    path = validation_path(data_dir, "Vedette")
    assert path.name == "vedette.tsv"
    key = f"chr2:{VEDETTE_START}-{VEDETTE_DEL[1]}:del"
    row = dict.fromkeys(VEDETTE_VALIDATION_HEADER, "")
    row.update(FID="1001", Variant=key, Sample="S3", Validation="present",
               Ignore="0", SiteId="vedette_site_1", CallerVersion=VEDETTE_CALLER)  # fmt: skip
    append_tsv_rows(path, VEDETTE_VALIDATION_HEADER, [list(row.values())])
    vmap = load_validation_map(path, "1001")
    assert vmap[(key, "S3")][0][0] == "present"
    merged = merge_validation_maps(vmap, {(key, "S3"): [("absent",)]})
    assert [x[0] for x in merged[(key, "S3")]] == ["present", "absent"]


# ---------------------------------------------------------------------------
# Deletions inside reads
# ---------------------------------------------------------------------------


def _small_del_reads(data_dir):
    _, windows = read_view_windows("del", "chr2", *VEDETTE_DEL, precise=True)
    cram = sample_dir(data_dir, "S3") / "sequences" / "S3.GRCh38_GIABv3.cram"
    [cells] = load_read_evidence([cram], windows, cram_reference(data_dir.parent))
    return windows, cells


def test_deletions_inside_reads_are_collected(data_dir):
    windows, [cell] = _small_del_reads(data_dir)
    assert len(windows) == 1  # 120 bp: one window
    [d] = cell.indels
    length = VEDETTE_DEL[1] - VEDETTE_DEL[0]
    assert (d.kind, d.x1, d.x2, d.length) == ("del", *VEDETTE_DEL, length)
    assert d.count == CRAM_INREAD_DELS
    # No pair, split or clip shows it.
    assert not cell.pairs and not cell.splits and not cell.clips


def test_deletions_inside_reads_give_the_breakpoints(data_dir):
    windows, cells = _small_del_reads(data_dir)
    start_call, end_call = suggest_breakpoints(
        cells, [None], [None], "del", *VEDETTE_DEL
    )
    # The reads give the last base kept (the VCF POS), one before the export's
    # start (the first base deleted), and the last base deleted.
    assert start_call.pos == VEDETTE_START - 1 and end_call.pos == VEDETTE_DEL[1]
    assert start_call.splits == CRAM_INREAD_DELS
    # Not for a duplication: an in-read deletion does not speak for one.
    dup_start, _ = suggest_breakpoints(cells, [None], [None], "dup", *VEDETTE_DEL)
    assert dup_start.pos is None


def test_deletions_inside_reads_are_drawn(data_dir):
    windows, cells = _small_del_reads(data_dir)
    option = read_option(
        [[(cells[0], None)]], ["S3 (carrier)"], windows=windows, breakpoints=[]
    )
    [series] = [s for s in option["series"] if s["name"] == SERIES_NAMES["in_read_del"]]
    [item] = series["data"]
    assert item["coords"][0][0] == VEDETTE_DEL[0]
    assert item["coords"][1][0] == VEDETTE_DEL[1]
    assert f"{CRAM_INREAD_DELS} reads with a 120 bp deletion" in item["tip"]
    assert SERIES_NAMES["in_read_del"] in option["legend"]["data"]
