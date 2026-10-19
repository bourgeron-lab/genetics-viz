"""Read-evidence (samplot-like) panel of the SV dialog.

Discordant pairs at their insert size, split reads, soft-clip clusters and the
copy number from the normalised coverage, in one plot per sample and window.
The windows are those of the IGV read view (``utils/sv_view.py``). CRAMs are
decoded off-thread with pysam; the copy number comes from the coverage bigWigs
and the cohort's control panel, like the Copy Number View.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

from nicegui import ui

from genetics_viz.components.offthread_panel import OffThreadPanel
from genetics_viz.components.sv_cn_panel import panel_search_dirs
from genetics_viz.config_model import get_reference_fasta
from genetics_viz.utils.cnv_coverage import CNRequest, load_cn_view
from genetics_viz.utils.sharding import get_sample_path
from genetics_viz.utils.sv_read_plots import read_chart_height, read_option
from genetics_viz.utils.sv_reads import load_read_evidence
from genetics_viz.utils.sv_view import Window, visibility_window_bp


def sample_cram(data_dir: Path, sample_id: str) -> Optional[Path]:
    path = (
        get_sample_path(data_dir, sample_id)
        / "sequences"
        / f"{sample_id}.GRCh38_GIABv3.cram"
    )
    return path if path.exists() else None


def load_read_view(
    data_dir: Path,
    cohort_dir: Optional[Path],
    cohort_name: str,
    params: Dict[str, Any],
    requests: Sequence[CNRequest],
    windows: Sequence[Window],
    reference: Optional[Path],
):
    """Reads (when a reference and narrow windows allow) and CN, per sample and window.

    Blocking: call it off-thread.
    """
    crams = [sample_cram(data_dir, r.sample_id) for r in requests]
    read_reads = reference is not None and all(
        w.width <= visibility_window_bp() for w in windows
    )
    reads = (
        load_read_evidence(crams, windows, reference)
        if read_reads
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
    rows = [
        [
            (reads[i][j] if reads[i] is not None else None, per_window[j][i])
            for j in range(len(windows))
        ]
        for i in range(len(requests))
    ]
    return panel, rows, read_reads


class ReadEvidencePanel:
    def __init__(
        self,
        *,
        data_dir: Path,
        cohort_dir: Optional[Path],
        cohort_name: str,
        params: Dict[str, Any],
        breakpoints: Sequence[int],
        label_for: Callable[[str], str],
    ) -> None:
        self.data_dir = data_dir
        self.cohort_dir = cohort_dir
        self.cohort_name = cohort_name
        self.params = params
        self.breakpoints = list(breakpoints)
        self.label_for = label_for
        self._panel = OffThreadPanel("read-panel")

    def refresh(self, requests: Sequence[CNRequest], windows: Sequence[Window]) -> None:
        requests, windows = list(requests), list(windows)
        reference = get_reference_fasta()
        self._panel.load(
            lambda: load_read_view(
                self.data_dir,
                self.cohort_dir,
                self.cohort_name,
                self.params,
                requests,
                windows,
                reference,
            ),
            lambda result: self._render(*result, requests, windows, reference),
            "Reading alignments…",
        )

    def _render(
        self,
        panel,
        rows,
        read_reads: bool,
        requests: List[CNRequest],
        windows: List[Window],
        reference: Optional[Path],
    ) -> None:
        if reference is None:
            ui.label(
                "Set 'reference_fasta' in the application config (the FASTA the "
                "CRAMs were aligned to) to see the reads."
            ).classes("text-sm text-amber-700 px-2")
        elif not read_reads:
            limit = visibility_window_bp() // 1000
            ui.label(
                f"A window is wider than {limit} kb: only the copy number is "
                "drawn. Switch to Split to see the reads at each breakpoint."
            ).classes("text-sm text-amber-700 px-2")
        has_data = [
            any(reads is not None or cn is not None for reads, cn in cells)
            for cells in rows
        ]
        shown = [(cells, r) for cells, r, ok in zip(rows, requests, has_data) if ok]
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
            breakpoints=self.breakpoints,
            panel=panel,
        )
        ui.echart(option).classes("w-full").style(
            f"height: {read_chart_height(len(shown))}px"
        ).mark("read-chart")
