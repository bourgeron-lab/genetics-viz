"""Read-evidence (samplot-like) panel of the SV dialog.

Pairs at their insert size, coloured by orientation, with split reads and a
coverage strip, for the same windows as the IGV read view
(``utils/sv_view.py``). The CRAMs are decoded off-thread with pysam.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, List, Optional, Sequence

from nicegui import ui

from genetics_viz.components.offthread_panel import OffThreadPanel
from genetics_viz.config_model import get_reference_fasta
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


class ReadEvidencePanel:
    def __init__(
        self,
        *,
        data_dir: Path,
        breakpoints: Sequence[int],
        label_for: Callable[[str], str],
    ) -> None:
        self.data_dir = data_dir
        self.breakpoints = list(breakpoints)
        self.label_for = label_for
        self._panel = OffThreadPanel("read-panel")

    def refresh(self, sample_ids: Sequence[str], windows: Sequence[Window]) -> None:
        reference = get_reference_fasta()
        if reference is None:
            self._panel.show_message(
                "Set 'reference_fasta' in the application config (the FASTA the "
                "CRAMs were aligned to) to see the read evidence."
            )
            return
        limit = visibility_window_bp()
        if any(w.width > limit for w in windows):
            self._panel.show_message(
                f"A window is wider than {limit // 1000} kb: switch to Split to see "
                "the reads at each breakpoint."
            )
            return
        sample_ids = list(sample_ids)
        windows = list(windows)
        crams = [sample_cram(self.data_dir, s) for s in sample_ids]
        self._panel.load(
            lambda: load_read_evidence(crams, windows, reference),
            lambda result: self._render(result, sample_ids),
            "Reading alignments…",
        )

    def _render(self, result, sample_ids: List[str]) -> None:
        shown = [(r, s) for r, s in zip(result, sample_ids) if r is not None]
        missing = [s for r, s in zip(result, sample_ids) if r is None]
        if missing:
            ui.label(f"No CRAM for: {', '.join(missing)}").classes(
                "text-sm text-gray-500 px-2"
            )
        if not shown:
            return
        option = read_option(
            [r for r, _ in shown],
            [f"{s} {self.label_for(s)}".strip() for _, s in shown],
            breakpoints=self.breakpoints,
        )
        ui.echart(option).classes("w-full").style(
            f"height: {read_chart_height(len(shown))}px"
        ).mark("read-chart")
