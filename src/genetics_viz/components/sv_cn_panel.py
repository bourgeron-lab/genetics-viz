"""Copy-number panel of the SV dialog: normalised coverage ÷ control panel.

Replaces the IGV by1000 bedgraph view when the carrier has a coverage bigWig
(``utils/cnv_coverage.py``). The bigWigs are read off-thread; the chart is
built once they are in, and rebuilt whenever the sample list changes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from nicegui import ui

from genetics_viz.components.offthread_panel import OffThreadPanel
from genetics_viz.utils.cnv_coverage import (
    CNRequest,
    cn_window,
    load_cn_view,
)
from genetics_viz.utils.cnv_plots import cn_chart_height, cn_option
from genetics_viz.utils.data import get_all_data_stores


def panel_search_dirs(data_dir: Path) -> List[Path]:
    """Where a panel path is looked for: the cohort's data directory, then the rest."""
    others = [Path(p) for p in get_all_data_stores() if Path(p) != data_dir]
    return [data_dir, *others]


class CopyNumberPanel:
    """A refreshable copy-number chart for one SV and a list of samples."""

    def __init__(
        self,
        *,
        data_dir: Path,
        cohort_dir: Optional[Path],
        cohort_name: str,
        params: Dict[str, Any],
        chrom: str,
        roi: Tuple[int, int],
        call: str,
        label_for: Callable[[str], str],
    ) -> None:
        self.data_dir = data_dir
        self.cohort_dir = cohort_dir
        self.cohort_name = cohort_name
        self.params = params
        self.chrom = chrom
        self.roi = roi
        self.call = call
        self.label_for = label_for
        self.window = cn_window(*roi)
        self._panel = OffThreadPanel("cn-panel")

    def _data_dirs(self) -> List[Path]:
        return panel_search_dirs(self.data_dir)

    def set_roi(self, roi: Tuple[int, int]) -> None:
        """New SV bounds (a curated position): the window follows them."""
        self.roi = roi
        self.window = cn_window(*roi)

    def refresh(self, requests: Sequence[CNRequest]) -> None:
        requests = list(requests)
        self._panel.load(
            lambda: load_cn_view(
                self.data_dir,
                self._data_dirs(),
                self.cohort_dir,
                self.cohort_name,
                self.params,
                requests,
                self.chrom,
                *self.window,
            ),
            lambda result: self._render(*result, requests),
            "Loading coverage…",
        )

    def _render(self, panel, tracks, requests: Sequence[CNRequest]) -> None:
        shown = [(t, r) for t, r in zip(tracks, requests) if t is not None]
        missing = [r.sample_id for t, r in zip(tracks, requests) if t is None]
        if panel is None:
            ui.label(
                "No CNV control panel for this cohort (cnv_panel_dir in its "
                "params file): showing each sample's own normalisation."
            ).classes("text-sm text-amber-700 px-2")
        if missing:
            ui.label(f"No normalised coverage for: {', '.join(missing)}").classes(
                "text-sm text-gray-500 px-2"
            )
        if not shown:
            return
        option = cn_option(
            [t for t, _ in shown],
            [self.label_for(r.sample_id) for _, r in shown],
            window=self.window,
            roi=self.roi,
            call=self.call,
            panel=panel,
        )
        ui.echart(option).classes("w-full").style(
            f"height: {cn_chart_height(len(shown))}px"
        ).mark("cn-chart")
