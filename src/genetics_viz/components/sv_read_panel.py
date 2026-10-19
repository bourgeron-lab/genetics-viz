"""Read-evidence (samplot-like) panel of the SV dialog.

One plot per sample and window shows:
- discordant pairs at their insert size;
- split reads and soft-clip clusters;
- the copy number from the normalised coverage (÷ panel);
- a finer read-depth line from the CRAM.

The windows are those of the IGV read view (``utils/sv_view.py``). The
carrier's clips and split reads suggest exact breakpoints, which can be taken
as the curated position.

CRAMs are decoded off-thread with pysam.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

from nicegui import ui

from genetics_viz.components.offthread_panel import OffThreadPanel
from genetics_viz.components.sv_cn_panel import panel_search_dirs
from genetics_viz.config_model import get_reference_fasta
from genetics_viz.utils.cnv_coverage import (
    CNRequest,
    autosomal_median_depth,
    find_sample_coverage,
    load_cn_view,
    pick_coverage,
)
from genetics_viz.utils.sharding import get_sample_path
from genetics_viz.utils.sv_breakpoints import SideSuggestion, suggest_breakpoints
from genetics_viz.utils.sv_read_plots import (
    Cell,
    Suggested,
    read_chart_height,
    read_option,
)
from genetics_viz.utils.sv_reads import load_read_evidence
from genetics_viz.utils.sv_view import Window, visibility_window_bp


def sample_cram(data_dir: Path, sample_id: str) -> Optional[Path]:
    path = (
        get_sample_path(data_dir, sample_id)
        / "sequences"
        / f"{sample_id}.GRCh38_GIABv3.cram"
    )
    return path if path.exists() else None


@dataclass
class ReadView:
    panel: Any
    rows: List[List[Cell]]
    reads_loaded: bool
    start_call: Optional[SideSuggestion]
    end_call: Optional[SideSuggestion]


def load_read_view(
    data_dir: Path,
    cohort_dir: Optional[Path],
    cohort_name: str,
    params: Dict[str, Any],
    requests: Sequence[CNRequest],
    windows: Sequence[Window],
    bounds: Sequence[int],
    reference: Optional[Path],
    include_low_mapq: bool = False,
    sv_type: str = "",
) -> ReadView:
    """Reads (when a reference and narrow windows allow), CN and read depth.

    Blocking: call it off-thread.
    """
    crams = [sample_cram(data_dir, r.sample_id) for r in requests]
    reads_loaded = reference is not None and all(
        w.width <= visibility_window_bp() for w in windows
    )
    reads = (
        load_read_evidence(crams, windows, reference, include_low_mapq=include_low_mapq)
        if reads_loaded
        else [None] * len(requests)
    )
    panel = None
    per_window = []
    for w in windows:
        panel, tracks = load_cn_view(
            data_dir,
            panel_search_dirs(data_dir),
            cohort_dir,
            cohort_name,
            params,
            requests,
            w.chrom,
            w.start,
            w.end,
        )
        per_window.append(tracks)

    medians = []
    for r in requests:
        cov = pick_coverage(find_sample_coverage(data_dir, r.sample_id), None)
        medians.append(autosomal_median_depth(cov) if cov else None)

    rows: List[List[Cell]] = []
    for i in range(len(requests)):
        row = []
        for j in range(len(windows)):
            cell_reads = reads[i][j] if reads[i] is not None else None
            depth = None
            if cell_reads is not None and medians[i]:
                scale = 2 / medians[i]
                depth = [[p, round(d * scale, 3)] for p, d in cell_reads.depth]
            row.append(Cell(cell_reads, per_window[j][i], depth))
        rows.append(row)

    start_call = end_call = None
    if rows:  # the carrier: reads when loaded, coverage either way
        carrier = rows[0]
        start_call, end_call = suggest_breakpoints(
            [c.reads for c in carrier],
            [c.cn for c in carrier],
            [c.depth for c in carrier],
            sv_type,
            *bounds,
        )
    return ReadView(panel, rows, reads_loaded, start_call, end_call)


_COVERAGE_WORDS = {
    "supported": "coverage agrees",
    "contradicted": "coverage disagrees",
    "unknown": "coverage unknown",
}


def _describe(call: Optional[SideSuggestion]) -> str:
    """One side in words: exact position, or the approximate one, and why."""
    if call is None:
        return "none"
    if call.pos is not None:
        parts = []
        if call.clips:
            parts.append(f"{call.clips} clipped")
        if call.splits:
            parts.append(f"{call.splits} split")
        if call.mismatched:
            parts.append(f"{call.mismatched} mis-oriented")
        text = f"{call.pos:,} {call.confidence} ({', '.join(parts)})"
        notes = []
        if call.coverage in _COVERAGE_WORDS:
            delta = (
                f" {call.coverage_delta:+.1f}"
                if call.coverage_delta is not None
                else ""
            )
            notes.append(_COVERAGE_WORDS[call.coverage] + delta)
        if call.outside_interval:
            notes.append("outside the pairs' interval")
        return text + (f", {', '.join(notes)}" if notes else "")
    if call.step is not None:
        text = (
            f"≈ {call.step.pos:,} (coverage step {call.step.delta:+.1f},"
            f" ±{call.step.resolution} bp)"
        )
    else:
        text = "none from reads or coverage"
    if call.interval is not None:
        lo, hi = call.interval
        text += f"; pairs place it in {lo:,}-{hi:,} ({call.interval_pairs} reads)"
    return text


def _plot_marks(view: ReadView) -> tuple[list[Suggested], list[tuple[int, int]]]:
    marks: list[Suggested] = []
    intervals: list[tuple[int, int]] = []
    for call in (view.start_call, view.end_call):
        if call is None:
            continue
        pos, approximate = call.position()
        if pos is not None:
            label = f"≈{pos:,}" if approximate else f"{pos:,} {call.confidence}"
            marks.append(Suggested(pos, approximate, label))
        if call.interval is not None:
            intervals.append(call.interval)
    return marks, intervals


class ReadEvidencePanel:
    def __init__(
        self,
        *,
        data_dir: Path,
        cohort_dir: Optional[Path],
        cohort_name: str,
        params: Dict[str, Any],
        label_for: Callable[[str], str],
        on_use_breakpoints: Callable[[Optional[int], Optional[int], bool], None],
    ) -> None:
        self.data_dir = data_dir
        self.cohort_dir = cohort_dir
        self.cohort_name = cohort_name
        self.params = params
        self.label_for = label_for
        self.on_use_breakpoints = on_use_breakpoints
        self._panel = OffThreadPanel("read-panel")

    def refresh(
        self,
        requests: Sequence[CNRequest],
        windows: Sequence[Window],
        bounds: Sequence[int],
        include_low_mapq: bool = False,
        sv_type: str = "",
    ) -> None:
        requests, windows, bounds = list(requests), list(windows), list(bounds)
        reference = get_reference_fasta()
        self._panel.load(
            lambda: load_read_view(
                self.data_dir,
                self.cohort_dir,
                self.cohort_name,
                self.params,
                requests,
                windows,
                bounds,
                reference,
                include_low_mapq,
                sv_type,
            ),
            lambda view: self._render(view, requests, windows, bounds, reference),
            "Reading alignments…",
        )

    def _render(
        self,
        view: ReadView,
        requests: List[CNRequest],
        windows: List[Window],
        bounds: List[int],
        reference: Optional[Path],
    ) -> None:
        if reference is None:
            ui.label(
                "Set 'reference_fasta' in the application config (the FASTA the "
                "CRAMs were aligned to) to see the reads."
            ).classes("text-sm text-amber-700 px-2")
        elif not view.reads_loaded:
            limit = visibility_window_bp() // 1000
            ui.label(
                f"A window is wider than {limit} kb: only the copy number is "
                "drawn. Switch to Split to see the reads at each breakpoint."
            ).classes("text-sm text-amber-700 px-2")

        marks, intervals = _plot_marks(view)
        start_call, end_call = view.start_call, view.end_call
        if start_call is not None and end_call is not None:
            (start_pos, start_approx) = start_call.position()
            (end_pos, end_approx) = end_call.position()
            with ui.column().classes("gap-0 px-2").mark("read-breakpoints"):
                for side, call in (("Start", start_call), ("End", end_call)):
                    with ui.row().classes("items-center gap-2"):
                        ui.icon("my_location", color="green").classes("text-base")
                        ui.label(f"{side}: {_describe(call)}").classes(
                            "text-sm text-green-800"
                            if call.pos is not None
                            else "text-sm text-gray-600"
                        )
                if start_pos is not None or end_pos is not None:
                    approximate = start_approx or end_approx
                    ui.button(
                        "Use approximate position"
                        if approximate
                        else "Use as curated position",
                        icon="content_paste",
                        on_click=lambda: self.on_use_breakpoints(
                            start_pos, end_pos, approximate
                        ),
                    ).props("dense outline color=green size=sm no-caps").classes(
                        "self-start mt-1"
                    )

        has_data = [
            any(c.reads is not None or c.cn is not None for c in cells)
            for cells in view.rows
        ]
        shown = [
            (cells, r) for cells, r, ok in zip(view.rows, requests, has_data) if ok
        ]
        missing = [r.sample_id for r, ok in zip(requests, has_data) if not ok]
        if missing:
            ui.label(f"No CRAM nor coverage for: {', '.join(missing)}").classes(
                "text-sm text-gray-500 px-2"
            )
        if not shown:
            return
        option = read_option(
            [cells for cells, _ in shown],
            [f"{r.sample_id} {self.label_for(r.sample_id)}".strip() for _, r in shown],
            windows=windows,
            breakpoints=bounds,
            suggested=marks,
            intervals=intervals,
            panel=view.panel,
        )
        ui.echart(option).classes("w-full").style(
            f"height: {read_chart_height(len(shown))}px"
        ).mark("read-chart")
