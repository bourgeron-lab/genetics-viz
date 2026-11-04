"""Vedette calls: parsing, carrier rows, filters and presets, link, own call.

Also the read evidence Vedette's small calls need: deletions inside reads.
"""

from datagen import (
    CRAM_INREAD_DELS,
    VEDETTE_DEL,
    VEDETTE_SOURCE,
    cram_reference,
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


def _rows(data_dir):
    df = v.read_vedette_table(v.family_vedette_file(data_dir, "1001"))
    return v.carrier_rows(df, labels={"S1": "father", "S2": "mother"})


# ---------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------


def test_impact_is_split_into_score_level_and_reason():
    assert v.parse_impact("+0.90 — élevé · DEL exonique RERE (NDD, s_het 0.1)") == (
        0.9,
        "high",
        "DEL exonique RERE (NDD, s_het 0.1)",
    )
    assert v.parse_impact("-0.60 — aucun gène codant") == (-0.6, "none", "")
    assert v.parse_impact("+0.00 — faible · intronique")[1] == "low"
    # A level Vedette adds later is kept as written, coloured grey.
    assert v.parse_impact("+0.20 — nouveau · x")[1] == "nouveau"
    assert v.level_info("nouveau")["color"] == v.UNKNOWN_COLOR
    assert v.parse_impact(None) == (None, "", "")


def test_inheritance_and_genes_are_parsed():
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


def test_one_row_per_carrier(data_dir):
    rows = _rows(data_dir)
    assert [(r["site_id"], r["sample"]) for r in rows] == [
        ("vedette_site_1", "S3"),
        ("vedette_site_2", "S1"),
        ("vedette_site_2", "S3"),
        ("vedette_site_3", "S1"),
        ("vedette_site_3", "S2"),
        ("vedette_site_3", "S3"),
        ("vedette_site_4", "S2"),
        ("vedette_site_4", "S3"),
    ]
    first = rows[0]
    assert first["_variant_key"] == f"chr2:{VEDETTE_DEL[0]}-{VEDETTE_DEL[1]}:del"
    assert (first["GT"], first["CN"], first["GQ"]) == ("0/1", 1, 45)
    assert first["inheritance"] == "de novo" and first["impact_level"] == "High"
    # Every member's genotype, the father named in the tooltip.
    assert [b["label"] for b in first["family_badges"]] == ["0/0·2", "0/0·2", "0/1·1"]
    assert first["family_badges"][0]["tooltip"].startswith("S1 (father): GT 0/0")
    # A founder carrier has no inheritance of its own.
    assert rows[1]["inheritance"] == "" and rows[1]["_inheritance"] == ""
    assert rows[-1]["_genes"] == ["XRG", "PLAIN", "SHANK3"]
    assert rows[-1]["genes"].endswith("(+2)")


def test_missing_and_unknown_columns_degrade(tmp_path):
    path = tmp_path / "F.vedette.tsv"
    path.write_text(
        "site_id\tcontig\tstart\tend\tsvtype\tnew_col\tGT_A\tGT_B\n"
        "s1\tchr1\t100\t300\tDEL\tx\t0/1\t0/0\n"
    )
    [row] = v.carrier_rows(v.read_vedette_table(path))
    assert row["sample"] == "A" and row["size"] == 200
    assert row["new_col"] == "x" and row["impact_score"] is None
    assert row["CN"] is None and row["cohort_af"] is None


def test_caller_version_is_the_last_source_line(data_dir):
    assert v.caller_version(data_dir, "1001") == VEDETTE_SOURCE
    assert v.caller_version(data_dir, "1002") == ""


# ---------------------------------------------------------------------------
# Filters, presets, link
# ---------------------------------------------------------------------------


def test_presets(data_dir):
    rows = _rows(data_dir)

    def keys(f, **kw):
        return [(r["site_id"], r["sample"]) for r in v.apply_filters(rows, f, **kw)]

    assert v.default_preset() == "rare"
    assert keys(v.VedetteFilters.preset("rare")) == [
        ("vedette_site_1", "S3"),
        ("vedette_site_4", "S2"),
        ("vedette_site_4", "S3"),
    ]
    assert keys(v.VedetteFilters.preset("denovo")) == [("vedette_site_1", "S3")]
    assert len(keys(v.VedetteFilters.preset("impact"))) == 3
    assert len(keys(v.VedetteFilters.preset("all"))) == 8
    # The page's member selection applies on top.
    assert keys(v.VedetteFilters.preset("all"), members=["S1"]) == [
        ("vedette_site_2", "S1"),
        ("vedette_site_3", "S1"),
    ]


def test_each_filter(data_dir):
    rows = _rows(data_dir)
    f = v.VedetteFilters.preset("all")
    f.types = ["DUP"]
    assert {r["site_id"] for r in v.apply_filters(rows, f)} == {"vedette_site_2"}
    f = v.VedetteFilters.preset("all")
    f.min_size, f.max_size = 100, 350
    assert {r["site_id"] for r in v.apply_filters(rows, f)} == {
        "vedette_site_1",
        "vedette_site_4",
    }
    f = v.VedetteFilters.preset("all")
    f.inheritance = [""]  # "not given": the founders' rows
    assert {r["sample"] for r in v.apply_filters(rows, f)} == {"S1", "S2"}
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
    custom.max_af, custom.min_size, custom.inheritance = 0.05, 200, ["maternal", ""]
    custom.genesets, custom.validations = ["tier:1"], ["present"]
    params, back = through(custom)
    assert params["vpre"] == "custom"
    assert (back.max_af, back.min_size, back.max_size) == (0.05, 200, None)
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
    own = v.find_own_call(vcf, "chr2", *VEDETTE_DEL, "DEL")
    assert own.call_id == "vedette_DEL_7" and own.passed and own.precise
    assert own.model == "v5n" and own.fields["CG"] == 4
    depth = v.find_own_call(vcf, "chrX", 2_800_000, 2_800_300, "DEL")
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


def test_vedette_validations_are_a_kind_of_their_own(data_dir):
    assert "Vedette" in VALIDATION_TYPES
    path = validation_path(data_dir, "Vedette")
    assert path.name == "vedette.tsv"
    key = f"chr2:{VEDETTE_DEL[0]}-{VEDETTE_DEL[1]}:del"
    row = dict.fromkeys(VEDETTE_VALIDATION_HEADER, "")
    row.update(FID="1001", Variant=key, Sample="S3", Validation="present",
               Ignore="0", SiteId="vedette_site_1", CallerVersion=VEDETTE_SOURCE)  # fmt: skip
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
    assert start_call.pos == VEDETTE_DEL[0] and end_call.pos == VEDETTE_DEL[1]
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
