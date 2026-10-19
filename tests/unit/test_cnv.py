"""CNV coverage, panel division, read windows and read evidence.

``utils/cnv_coverage.py``, ``utils/sv_view.py``, ``utils/sv_reads.py`` and
their option builders, on the synthetic trio of ``datagen``.
"""

import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from datagen import (
    CNV_LOSS,
    CNV_MASKED,
    CNV_PANEL_CLUSTER_PATH,
    CNV_TAG,
    CNV_MEDIAN_DEPTH,
    CRAM_CLIPPED,
    CRAM_CONCORDANT_LEFT,
    CRAM_SPLITS,
    PAR1_END,
    cram_reference,
    sample_dir,
)
from genetics_viz.utils import cnv_coverage as cnv
from genetics_viz.utils.cnv_plots import cn_option
from genetics_viz.utils.pipeline_params import load_params
from genetics_viz.utils.sv_breakpoints import coverage_profile, suggest_breakpoints
from genetics_viz.utils.sv_read_plots import (
    SERIES_NAMES,
    Cell,
    Suggested,
    _arc,
    read_option,
)
from genetics_viz.utils.sv_reads import load_read_evidence
from genetics_viz.utils.sv_view import auto_mode, read_view_windows


def _panel(data_dir: Path) -> cnv.Panel:
    params = load_params(data_dir / "cohorts" / "DEMO" / "DEMO.params.yml")
    panel = cnv.resolve_cohort_panel(params, [data_dir])
    assert panel is not None
    return panel


def _values(track: cnv.CNTrack, lo: int, hi: int, *, masked=False):
    return [
        v
        for p, v, m in zip(track.positions, track.values, track.masked)
        if lo <= p < hi and m == masked and v is not None
    ]


# ---------------------------------------------------------------------------
# Files and panel
# ---------------------------------------------------------------------------


def test_sample_coverage_is_found_by_tag(data_dir):
    found = cnv.find_sample_coverage(data_dir, "S3")
    assert [c.tag for c in found] == [CNV_TAG]
    assert cnv.find_sample_coverage(data_dir, "S6") == []
    assert cnv.pick_coverage(found, "cov_b1000q0") is found[0]  # falls back


def test_panel_dir_resolves_the_cluster_path_tail(data_dir, tmp_path):
    other = tmp_path / "other-data"
    other.mkdir()
    resolved = cnv.panel_dir_from_params(CNV_PANEL_CLUSTER_PATH, [other, data_dir])
    assert resolved == data_dir / "cohorts" / "DEMO" / "cnv_panel" / "v1"
    assert cnv.panel_dir_from_params("/nowhere/cnv_panel/v1", [data_dir]) is None
    assert cnv.panel_dir_from_params("", [data_dir]) is None
    assert cnv.panel_dir_from_params(None, [data_dir]) is None


def test_a_controls_cohort_uses_the_panel_it_built(data_dir):
    demo = data_dir / "cohorts" / "DEMO"
    # Its params name no cnv_panel_dir, only the version it builds.
    for params in ({"cnv_panel_version": "1"}, {}):
        panel = cnv.resolve_cohort_panel(params, [data_dir], demo)
        assert panel is not None and panel.version == "1"
    assert cnv.resolve_cohort_panel({}, [data_dir]) is None  # no cohort dir given


def test_panel_json_is_read(data_dir):
    panel = _panel(data_dir)
    assert (panel.version, panel.tag, panel.bin) == ("1", CNV_TAG, 250)
    assert panel.label == "DEMO v1"
    # 1-based inclusive in the file, half-open 0-based here.
    assert panel.par["chrX"][0] == (10000, PAR1_END)


def test_a_panel_without_its_median_is_ignored(data_dir):
    panel = _panel(data_dir)
    panel.bigwig("autosomes", "median").unlink()
    (panel.directory / "panel.json").write_text(
        (panel.directory / "panel.json").read_text() + " "
    )  # new stat, so the cache does not answer
    assert cnv.load_panel(panel.directory) is None


def test_fit_rows_and_karyotype(data_dir):
    panel = _panel(data_dir)
    fits = cnv.read_cohort_fit(data_dir / "cohorts" / "DEMO", "DEMO", panel)
    assert fits[("1001", "S3")].fit_mad == pytest.approx(0.41)
    assert cnv.karyotype(fits[("1001", "S2")], None) == "XX"
    assert cnv.karyotype(None, "1.0") == "XY"
    assert cnv.karyotype(None, "-9") is None
    odd = cnv.FitRow("f", "s", "male", "XXY", None, None, None)
    assert cnv.karyotype(odd, "1") is None


# ---------------------------------------------------------------------------
# Copy number
# ---------------------------------------------------------------------------


def test_the_carriers_loss_reads_one_copy(data_dir):
    panel = _panel(data_dir)
    cov = cnv.pick_coverage(cnv.find_sample_coverage(data_dir, "S3"), panel.tag)
    track = cnv.copy_number(cov, panel, "XY", "chr2", 165_000_000, 165_500_000)
    assert track.mode == "corrected" and track.reason is None
    inside = _values(track, CNV_LOSS[0] + 1000, CNV_LOSS[1] - 1000)
    outside = _values(track, 165_000_000, CNV_LOSS[0] - 1000)
    assert inside and np.allclose(inside, 1.0)
    assert outside and np.allclose(outside, 2.0)


def test_unusable_panel_bins_are_masked(data_dir):
    panel = _panel(data_dir)
    cov = cnv.pick_coverage(cnv.find_sample_coverage(data_dir, "S1"), panel.tag)
    track = cnv.copy_number(cov, panel, "XY", "chr2", *CNV_MASKED)
    assert all(track.masked)
    assert _values(track, *CNV_MASKED, masked=True)  # still drawn, greyed


@pytest.mark.parametrize(
    ("sample", "kary", "par", "non_par"),
    [("S1", "XY", 2.0, 1.0), ("S2", "XX", 2.0, 2.0)],
)
def test_chrx_is_divided_by_the_sex_group(data_dir, sample, kary, par, non_par):
    panel = _panel(data_dir)
    cov = cnv.pick_coverage(cnv.find_sample_coverage(data_dir, sample), panel.tag)
    track = cnv.copy_number(cov, panel, kary, "chrX", 2_700_000, 2_900_000)
    assert track.mode == "corrected"
    assert np.allclose(_values(track, 2_700_000, PAR1_END - 500), par)
    assert np.allclose(_values(track, PAR1_END + 500, 2_900_000), non_par)


def test_unknown_sex_is_not_divided_on_x(data_dir):
    panel = _panel(data_dir)
    cov = cnv.pick_coverage(cnv.find_sample_coverage(data_dir, "S1"), panel.tag)
    track = cnv.copy_number(cov, panel, None, "chrX", 2_700_000, 2_900_000)
    assert track.mode == "raw" and "karyotype unknown" in track.reason


def test_no_panel_or_another_tag_shows_the_raw_value(data_dir):
    panel = _panel(data_dir)
    cov = cnv.pick_coverage(cnv.find_sample_coverage(data_dir, "S3"), panel.tag)
    raw = cnv.copy_number(cov, None, "XY", "chr2", 165_000_000, 165_500_000)
    assert raw.mode == "raw" and "no control panel" in raw.reason
    other = dataclasses.replace(panel, tag="cov_b1000q1")
    mismatch = cnv.copy_number(cov, other, "XY", "chr2", 165_000_000, 165_500_000)
    assert mismatch.mode == "raw" and "cov_b1000q1" in mismatch.reason
    assert np.allclose(_values(raw, CNV_LOSS[0] + 1000, CNV_LOSS[1] - 1000), 1.0)


def test_wide_regions_are_aggregated(data_dir):
    panel = _panel(data_dir)
    cov = cnv.pick_coverage(cnv.find_sample_coverage(data_dir, "S3"), panel.tag)
    track = cnv.copy_number(
        cov, panel, "XY", "chr2", 164_800_000, 165_700_000, cfg={"display_points": 100}
    )
    assert len(track.positions) <= 100 and track.bin_bp >= 9000


def test_load_cn_view_reads_a_trio(data_dir):
    params = load_params(data_dir / "cohorts" / "DEMO" / "DEMO.params.yml")
    reqs = [
        cnv.CNRequest("S3", "1001", "1"),
        cnv.CNRequest("S1", "1001", "1"),
        cnv.CNRequest("S6", "AU-X-1", "2"),
    ]
    panel, tracks = cnv.load_cn_view(
        data_dir,
        [data_dir],
        data_dir / "cohorts" / "DEMO",
        "DEMO",
        params,
        reqs,
        "chr2",
        *cnv.cn_window(*CNV_LOSS),
    )
    assert panel is not None
    assert tracks[2] is None  # S6 has no coverage file
    assert tracks[0].fit_mad == pytest.approx(0.41) and tracks[0].karyotype == "XY"

    option = cn_option(
        [t for t in tracks if t],
        ["(carrier)", "(father)"],
        window=cnv.cn_window(*CNV_LOSS),
        roi=CNV_LOSS,
        call="loss",
        panel=panel,
    )
    assert len(option["grid"]) == 2
    assert option["dataZoom"][0]["xAxisIndex"] == [0, 1]
    area = option["series"][0]["markArea"]["data"][0]
    assert (area[0]["xAxis"], area[1]["xAxis"]) == CNV_LOSS
    # S3's fit is worse than the panel's own cut-off: flagged.
    assert option["title"][0]["textStyle"]["color"] == "#b45309"


def test_cn_window_keeps_flanks_of_a_short_event():
    lo, hi = cnv.cn_window(1_000_000, 1_000_500)
    assert hi - lo == 20_000 and lo < 1_000_000 < 1_000_500 < hi


# ---------------------------------------------------------------------------
# Read windows
# ---------------------------------------------------------------------------


def test_short_svs_get_one_window_and_long_ones_two():
    assert auto_mode("del", 1_000_000, 1_010_000) == "single"
    assert auto_mode("del", 1_000_000, 1_100_000) == "split"
    assert auto_mode("ins", 1_000_000, 1_100_000) == "single"
    assert auto_mode("bnd", 1_000_000, 1_000_001) == "split"

    mode, (one,) = read_view_windows("del", "chr1", 1_000_000, 1_000_400)
    assert mode == "single" and one.width == 3000  # never under two flanks
    mode, (left, right) = read_view_windows("dup", "chr1", 1_000_000, 2_000_000)
    assert (left.start, left.end, right.start, right.end) == (
        998_500,
        1_001_500,
        1_998_500,
        2_001_500,
    )


def test_the_mode_can_be_forced_and_a_mate_placed():
    mode, windows = read_view_windows("del", "chr1", 1_000_000, 1_010_000, mode="split")
    assert mode == "split" and len(windows) == 2
    _, (_, mate) = read_view_windows("bnd", "chr1", 5000, 5001, mate=("chr7", 9000))
    assert mate.chrom == "chr7" and mate.start == 7500


# ---------------------------------------------------------------------------
# Read evidence
# ---------------------------------------------------------------------------


def _s3_reads(data_dir, cfg=None, include_low_mapq=False):
    _, windows = read_view_windows("del", "chr2", *CNV_LOSS)
    cram = sample_dir(data_dir, "S3") / "sequences" / "S3.GRCh38_GIABv3.cram"
    (cells, missing) = load_read_evidence(
        [cram, None],
        windows,
        cram_reference(data_dir.parent),
        cfg=cfg,
        include_low_mapq=include_low_mapq,
    )
    assert missing is None
    return windows, cells


def test_read_signatures_are_classed(data_dir):
    _, (left, right) = _s3_reads(data_dir)
    got = {}
    for p in left.pairs:
        got[p.cls] = got.get(p.cls, 0) + 1
    # Concordant pairs are only counted, for the expected insert band.
    assert got == {"deletion": 6, "inversion_pp": 2, "duplication": 1}
    # The fixture's pairs reach past the window; those wholly outside it are
    # not fetched.
    assert 0.9 * CRAM_CONCORDANT_LEFT < left.n_concordant <= CRAM_CONCORDANT_LEFT
    assert left.insert_lo < left.insert_median < left.insert_hi < left.insert_max
    assert [m.mate for m in left.interchrom] == ["chr5:1000001"]
    # The deletion pairs reach across to the right breakpoint, read bodies at both ends.
    for p in left.pairs:
        if p.cls == "deletion":
            assert p.x2 > CNV_LOSS[1] and p.x1 < p.left_end <= p.right_start < p.x2
    assert sum(p.cls == "deletion" for p in right.pairs) == 6


def test_split_reads_keep_both_parts(data_dir):
    _, (left, _) = _s3_reads(data_dir)
    assert len(left.splits) == CRAM_SPLITS
    split = left.splits[0]
    assert (split.a1, split.a2) == (CNV_LOSS[0] - 60, CNV_LOSS[0])
    assert (split.b1, split.b2) == (CNV_LOSS[1], CNV_LOSS[1] + 40)
    assert split.y == CNV_LOSS[1] - CNV_LOSS[0] and split.cls == "split"
    # The read is clipped after its part, the SA part before its own: a
    # deletion junction, from the end of one to the start of the other.
    assert (split.ja, split.ja_side, split.jb, split.jb_side) == (
        CNV_LOSS[0],
        "right",
        CNV_LOSS[1],
        "left",
    )


def test_soft_clips_cluster_at_the_breakpoint(data_dir):
    _, (left, _) = _s3_reads(data_dir)
    # Three clipped reads and the two split reads end at the breakpoint.
    assert [(c.pos, c.side, c.count) for c in left.clips] == [
        (CNV_LOSS[0], "right", CRAM_CLIPPED + CRAM_SPLITS)
    ]
    _, (strict, _) = _s3_reads(data_dir, cfg={"min_clip_reads": 10})
    assert strict.clips == []


def test_read_option_draws_one_grid_per_cell_with_two_axes(data_dir):
    windows, (left, right) = _s3_reads(data_dir)
    panel = _panel(data_dir)
    cov = cnv.pick_coverage(cnv.find_sample_coverage(data_dir, "S3"), panel.tag)
    cns = [cnv.copy_number(cov, panel, "XY", w.chrom, w.start, w.end) for w in windows]
    option = read_option(
        [[(left, cns[0]), (right, cns[1])], [(None, cns[0]), (None, None)]],
        ["S3 (carrier)", "S1 (father)"],
        windows=windows,
        breakpoints=list(CNV_LOSS),
        panel=panel,
    )
    assert len(option["grid"]) == 4  # 2 samples x 2 windows
    assert [a["type"] for a in option["yAxis"][:2]] == ["log", "value"]
    assert option["yAxis"][1]["position"] == "right"
    # Two rulers and two drag-only zooms, one per window; no wheel zoom.
    kinds = [z["type"] for z in option["dataZoom"]]
    assert kinds.count("slider") == 2 and kinds.count("inside") == 2
    assert all(
        not z["zoomOnMouseWheel"] for z in option["dataZoom"] if z["type"] == "inside"
    )
    by_name = {}
    for s in option["series"]:
        by_name.setdefault(s["name"], []).append(s)
    assert SERIES_NAMES["copy_number"] in by_name
    # The carrier's CN on the right axis reads 1 inside the loss at the right window.
    right_cn = [
        s
        for s in by_name[SERIES_NAMES["copy_number"]]
        if s["xAxisIndex"] == 1 and s["yAxisIndex"] == 3
    ][0]
    inside = [v for x, v in right_cn["data"] if v is not None and x < CNV_LOSS[1] - 300]
    assert inside and all(abs(v - 1.0) < 0.01 for v in inside)
    # Split reads: parts plus an arc rising above them.
    arc = next(s for s in by_name[SERIES_NAMES["split"]] if s.get("polyline"))
    coords = arc["data"][0]["coords"]
    assert max(y for _, y in coords) > coords[0][1]
    assert SERIES_NAMES["soft_clip"] in by_name
    band = by_name[SERIES_NAMES["insert_band"]][0]["markArea"]["data"][0]
    assert band[0]["yAxis"] == left.insert_lo
    assert "no CRAM" in option["title"][3]["text"]  # 2 window titles, then rows


def test_low_mapq_evidence_is_kept_only_when_asked(data_dir):
    _, (left, _) = _s3_reads(data_dir)
    assert not any(p.low for p in left.pairs)
    _, (with_low, _) = _s3_reads(data_dir, include_low_mapq=True)
    low = [p for p in with_low.pairs if p.low]
    assert [p.cls for p in low] == ["deletion"]
    # The low pair does not move the expected insert band.
    assert with_low.insert_lo == left.insert_lo


def test_read_depth_is_binned_from_the_reads(data_dir):
    _, (left, _) = _s3_reads(data_dir)
    assert len(left.depth) == left.window.width // 25
    assert max(d for _, d in left.depth) > 0


def _s3_evidence(data_dir, cn_sample="S3", strip_reads=False):
    """S3's read cells, CN tracks and scaled depth, as the panel builds them."""
    windows, cells = _s3_reads(data_dir)
    if strip_reads:
        cells = [dataclasses.replace(c, clips=[], splits=[]) for c in cells]
    panel = _panel(data_dir)
    cov = cnv.pick_coverage(cnv.find_sample_coverage(data_dir, cn_sample), panel.tag)
    cns = [cnv.copy_number(cov, panel, "XY", w.chrom, w.start, w.end) for w in windows]
    depths = (
        [[[p, d * 2 / CNV_MEDIAN_DEPTH] for p, d in c.depth] for c in cells]
        if cn_sample == "S3"
        else [None, None]
    )
    return cells, cns, depths


def test_the_fixture_depth_reads_as_copy_number(data_dir):
    cells, cns, depths = _s3_evidence(data_dir)
    profile = coverage_profile(depths[0], cns[0])
    outside = [v for p, v in profile if CNV_LOSS[0] - 1200 < p < CNV_LOSS[0] - 100]
    inside = [v for p, v in profile if CNV_LOSS[0] + 100 < p < CNV_LOSS[0] + 1200]
    assert 1.7 < sum(outside) / len(outside) < 2.3
    assert 0.8 < sum(inside) / len(inside) < 1.3


def test_a_deletion_is_suggested_from_clips_splits_and_coverage(data_dir):
    cells, cns, depths = _s3_evidence(data_dir)
    start_call, end_call = suggest_breakpoints(cells, cns, depths, "del", *CNV_LOSS)
    # Start: three clipped reads, the two split reads' clipped ends, and their
    # junction; coverage drops there, and the deletion pairs allow it.
    assert (start_call.pos, start_call.clips, start_call.splits) == (
        CNV_LOSS[0],
        CRAM_CLIPPED + CRAM_SPLITS,
        CRAM_SPLITS,
    )
    assert start_call.coverage == "supported" and start_call.coverage_delta < -0.5
    lo, hi = start_call.interval
    assert lo <= CNV_LOSS[0] <= hi and start_call.interval_pairs == 6
    assert start_call.confidence == "strong"
    # End: only the split reads' other part (2 x 2), coverage rises back.
    assert (end_call.pos, end_call.clips, end_call.splits) == (CNV_LOSS[1], 0, 2)
    assert end_call.coverage == "supported" and end_call.confidence == "moderate"
    assert end_call.position() == (CNV_LOSS[1], False)


def test_a_threshold_too_high_leaves_a_side_without_an_exact_position(data_dir):
    cells, cns, depths = _s3_evidence(data_dir)
    _, end_call = suggest_breakpoints(
        cells, cns, depths, "del", *CNV_LOSS, cfg={"min_breakpoint_score": 5}
    )
    assert end_call.pos is None and end_call.confidence is None


def test_mis_oriented_evidence_counts_less(data_dir):
    cells, cns, depths = _s3_evidence(data_dir)
    # Deletion-shaped clips and junctions do not make a duplication's start,
    # and the coverage does not rise there either.
    start_call, _ = suggest_breakpoints(cells, cns, depths, "dup", *CNV_LOSS)
    assert start_call.pos is None and start_call.step is None


def test_an_inversion_is_judged_on_reads_alone(data_dir):
    cells, cns, depths = _s3_evidence(data_dir)
    start_call, _ = suggest_breakpoints(cells, cns, depths, "inv", *CNV_LOSS)
    # Clips on either side are fine for an inversion; same-strand split reads
    # are not its kind, so they count for a quarter.
    assert start_call.pos == CNV_LOSS[0] and start_call.mismatched == CRAM_SPLITS
    assert start_call.coverage == "n/a" and start_call.step is None
    # The fixture's ++ pairs sit too far left for a junction there: weak.
    assert start_call.outside_interval and start_call.confidence == "weak"


def test_without_reads_the_coverage_step_is_approximate(data_dir):
    cells, cns, depths = _s3_evidence(data_dir, strip_reads=True)
    start_call, end_call = suggest_breakpoints(cells, cns, depths, "del", *CNV_LOSS)
    assert start_call.pos is None and start_call.coverage == "unknown"
    pos, approximate = start_call.position()
    # The fixture drops whole pairs on the deleted haplotype, which blurs its
    # depth step over about a fragment (~350 bp) on the loss side.
    assert approximate and CNV_LOSS[0] - 100 <= pos <= CNV_LOSS[0] + 400
    assert start_call.step.delta < -0.5
    # The pairs still place it.
    lo, hi = start_call.interval
    assert lo <= CNV_LOSS[0] <= hi


def test_a_flat_coverage_contradicts_a_deletion(data_dir):
    # S3's reads, but the father's flat CN and no read depth.
    cells, cns, depths = _s3_evidence(data_dir, cn_sample="S1")
    start_call, _ = suggest_breakpoints(cells, cns, depths, "del", *CNV_LOSS)
    assert start_call.pos == CNV_LOSS[0]
    assert start_call.coverage == "contradicted" and start_call.confidence == "weak"


def test_no_evidence_suggests_nothing():
    start_call, end_call = suggest_breakpoints(
        [None, None], [None, None], [None, None], "del", *CNV_LOSS
    )
    for call in (start_call, end_call):
        assert call.position() == (None, False) and call.interval is None


def test_split_read_arcs_are_half_ellipses():
    points = _arc(1000, 2000, 500, 4)
    xs, ys = zip(*points)
    assert ys[0] == ys[-1] == 500 and max(ys) == pytest.approx(2000)
    # Steep ends: the first step rises a lot over very little x.
    assert (ys[1] / ys[0] - 1) / ((xs[1] - xs[0]) / 1000) > 10


def test_read_option_draws_depth_and_suggested_breakpoints(data_dir):
    windows, (left, right) = _s3_reads(data_dir, include_low_mapq=True)
    depth = [[p, d * 2 / 30] for p, d in left.depth]
    option = read_option(
        [[Cell(left, None, depth), Cell(right, None, None)]],
        ["S3 (carrier)"],
        windows=windows,
        breakpoints=list(CNV_LOSS),
        suggested=[
            Suggested(CNV_LOSS[0] + 10, False, "exact"),
            Suggested(CNV_LOSS[1] + 10, True, "≈"),
        ],
        intervals=[(CNV_LOSS[0] - 50, CNV_LOSS[0] + 50)],
    )
    names = [s["name"] for s in option["series"]]
    assert SERIES_NAMES["read_depth"] in names
    marks = option["series"][0]["markLine"]["data"]
    green = [m for m in marks if "label" in m]
    assert [m["xAxis"] for m in green] == [CNV_LOSS[0] + 10]  # left window only
    assert green[0]["lineStyle"]["type"] == "solid"
    right_marks = [
        m
        for s in option["series"]
        if s.get("markLine") and s["xAxisIndex"] == 1
        for m in s["markLine"]["data"]
        if "label" in m
    ]
    assert right_marks[0]["lineStyle"]["type"] == "dashed"
    areas = option["series"][0]["markArea"]["data"]
    assert [a[0].get("xAxis") for a in areas if "xAxis" in a[0]] == [CNV_LOSS[0] - 50]
    faded = [
        s
        for s in option["series"]
        if s["name"] == SERIES_NAMES["deletion"] and s["lineStyle"]["opacity"] < 0.5
    ]
    assert faded  # the MAPQ-0 pair


def test_imprecise_calls_get_wider_breakpoint_windows():
    _, (left, right) = read_view_windows(
        "del", "chr1", 1_000_000, 2_000_000, precise=False
    )
    assert (left.start, left.end) == (995_000, 1_005_000)
    assert right.width == 10_000


def test_coverage_json_is_written_per_sample(data_dir):
    meta = json.loads(
        (sample_dir(data_dir, "S1") / "coverage" / f"S1.{CNV_TAG}.json").read_text()
    )
    assert meta["tag"] == CNV_TAG
    cov = cnv.find_sample_coverage(data_dir, "S1")[0]
    assert cnv.autosomal_median_depth(cov) == CNV_MEDIAN_DEPTH
