"""The cohort SVs tab: the lazy scan, filters, aneuploidies, recurrent regions."""

import datagen
from genetics_viz.utils import cohort_svs as cs

HEADER = ["chr", "start", "end", "ratio", "zscore", "type", "barcode", "exonic_symbol"]


def _bed(tmp_path, rows):
    return datagen.write_tsv(tmp_path / "C_aberrations.bed", HEADER, rows)


def test_scan_keeps_only_calls_past_the_permissive_thresholds(tmp_path):
    path = _bed(
        tmp_path,
        [
            ("2", 1_000_000, 3_000_000, -1.0, -9, "primary", "S1", "SCN2A"),
            ("chr2", 1_000_000, 3_000_000, 0.22, 2.6, "primary", "S2", ""),
            ("chr2", 1_000_000, 3_000_000, 0.05, 0.4, "primary", "S3", ""),  # below
            ("chr2", 1_000_000, 1_050_000, -1.0, -9, "primary", "S4", ""),  # small
        ],
    )
    calls = cs.scan_calls(path)
    assert calls["sample"].to_list() == ["S1", "S2"]
    assert calls["chr"].to_list() == ["chr2", "chr2"]
    assert calls["strength"].to_list() == ["robust", "permissive"]
    rows = cs.call_rows(calls, {"S1": ["F1"], "S2": ["F2", "F9"]})
    # S2 sits in two families: one row each.
    assert [(r["FID"], r["sample"]) for r in rows] == [
        ("F1", "S1"),
        ("F2", "S2"),
        ("F9", "S2"),
    ]
    assert rows[0]["Variant"] == "chr2:1000000-3000000:del"
    assert [r["sample"] for r in cs.filter_calls(rows, None)] == ["S1"]
    assert len(cs.filter_calls(rows, None, strengths=cs.STRENGTHS)) == 3
    assert cs.filter_calls(rows, {("F2", "S2")}, strengths=()) == [rows[1]]


def test_whole_chromosome_gain_is_an_aneuploidy(tmp_path):
    # chr21's q arm is about 12.0-46.7 Mb: a gain over all of it.
    path = _bed(tmp_path, [("chr21", 12_000_000, 46_700_000, 0.5, 8, "p", "S1", "")])
    rows = cs.call_rows(cs.scan_calls(path), {"S1": ["F1"]})
    events = cs.aneuploidy_rows(rows, 0.75)
    assert [(e["FID"], e["sample"]) for e in events] == [("F1", "S1")]
    assert events[0]["events"].startswith("+21")


def test_recurrent_regions_merge_overlapping_calls_of_two_families():
    def row(fid, sample, start, end, direction="loss", genes=()):
        return {"FID": fid, "sample": sample, "chr": "chr2", "start": start, "end": end,
                "direction": direction, "_genes": list(genes)}  # fmt: skip

    rows = [
        row("F1", "S1", 1_000, 5_000, genes=["A"]),
        row("F2", "S2", 4_000, 9_000, genes=["B"]),
        row("F2", "S2", 4_000, 9_000, direction="gain"),
        row("F3", "S3", 20_000, 30_000),
    ]
    regions = cs.recurrent_regions(rows)
    assert [(g["region"], g["n_families"], g["exonic_symbol"]) for g in regions] == [
        ("chr2:1000-9000", 2, "A,B")
    ]
