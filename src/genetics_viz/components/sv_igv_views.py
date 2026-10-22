"""SV dialog: the igv.js views (CNV coverage bedgraph, VAF, read level).

Each view is one igv.js browser in a div, rebuilt for the current samples.
``SVViewState`` is what they share with the rest of the dialog: the SV, the
samples shown and the current bounds, which the read views follow.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from nicegui import ui

from genetics_viz.utils.data import get_static_prefix
from genetics_viz.utils.sharding import get_sample_path, get_sample_url
from genetics_viz.utils.sv_view import Window, read_view_windows, visibility_window_bp


@dataclass
class SVViewState:
    """The SV and its samples, as every view of the dialog draws them."""

    data_dir: Path
    chrom: str
    #: The called position: titles, ROI names and validation keys use it.
    start: str
    end: str
    sample: str
    sv_type: str
    call: Any
    #: The locus and ROI the views open on: curated when a validation set one.
    igv_locus: str
    roi_start: str
    roi_end: str
    label_for: Callable[[str], str]
    #: Samples shown besides the carrier, in the order added.
    additional: List[str] = field(default_factory=list)
    #: The SV's current bounds: curated when set, and updated by every
    #: curated-position change. The read windows follow them.
    bounds: Dict[str, int] = field(default_factory=dict)
    #: One window around a short SV, one per breakpoint for a long one
    #: (utils/sv_view.py). Shared by the IGV reads and the read panel.
    read_mode: str = "auto"
    low_mapq: bool = False

    def __post_init__(self) -> None:
        if not self.bounds:
            self.bounds = {"start": int(self.roi_start), "end": int(self.roi_end)}

    def samples(self) -> List[str]:
        return [self.sample, *self.additional]

    def label(self, sample_id: str) -> str:
        return f"{sample_id} {self.label_for(sample_id)}".strip()

    @property
    def roi_color(self) -> str:
        return (
            "rgba(0,128,0,0.2)"  # Green for gain
            if "GAIN" in str(self.call).upper()
            else "rgba(255,0,0,0.2)"  # Red for loss
        )

    def read_windows(self) -> Tuple[str, List[Window]]:
        # WisecondorX breakpoints sit on 1 kb bins: wider windows.
        return read_view_windows(
            self.sv_type,
            self.chrom,
            self.bounds["start"],
            self.bounds["end"],
            mode=self.read_mode,
            precise=False,
        )

    def describe_read_mode(self, used_mode: str) -> str:
        return (
            f"{'auto: ' if self.read_mode == 'auto' else ''}"
            + ("one window" if used_mode == "single" else "one window per breakpoint")
            + f" (single up to {visibility_window_bp() // 1000} kb shown)"
        )


def _roi_js(state: SVViewState, name_span: str, start: Any, end: Any) -> str:
    """igv.js ``roi`` option: the SV's region, named by ``name_span``."""
    return f"""[
        {{
            name: "{state.chrom}:{name_span} - {state.call}",
            color: "{state.roi_color}",
            features: [
                {{
                    chr: "{state.chrom}",
                    start: {start},
                    end: {end},
                    name: "{state.call}"
                }}
            ]
        }}
    ]"""


def _sample_file(state: SVViewState, sample_id: str, name: str) -> Optional[str]:
    """The URL of ``sequences/<name>`` for a sample, if the file exists."""
    path = get_sample_path(state.data_dir, sample_id) / "sequences" / name
    if not path.exists():
        return None
    return f"{get_static_prefix()}/{get_sample_url(state.data_dir, sample_id)}/sequences/{name}"


def has_bedgraph(data_dir: Path, sample_id: str) -> bool:
    return (
        get_sample_path(data_dir, sample_id)
        / "sequences"
        / f"{sample_id}.by1000.bedgraph.gz"
    ).exists()


def coverage_tracks(state: SVViewState) -> List[Dict[str, Any]]:
    """The by1000 coverage bedgraph of each sample that has one."""
    tracks = []
    for sid in state.samples():
        url = _sample_file(state, sid, f"{sid}.by1000.bedgraph.gz")
        if url:
            tracks.append(
                {
                    "name": state.label(sid),
                    "type": "wig",
                    "format": "bedgraph",
                    "url": url,
                    "indexURL": f"{url}.tbi",
                    "height": 100,
                    "autoscaleGroup": "cnv",
                }
            )
    return tracks


def vaf_tracks(state: SVViewState) -> List[Dict[str, Any]]:
    """The VAF bedgraph of each sample that has one, as points in [0, 1]."""
    tracks = []
    for sid in state.samples():
        url = _sample_file(state, sid, f"{sid}.vaf.bedgraph.gz")
        if url:
            tracks.append(
                {
                    "name": state.label(sid),
                    "type": "wig",
                    "format": "bedgraph",
                    "graphType": "points",
                    "url": url,
                    "indexURL": f"{url}.tbi",
                    "height": 100,
                    "min": 0,
                    "max": 1,
                    "autoscale": False,
                }
            )
    return tracks


def cram_tracks(state: SVViewState) -> List[Dict[str, Any]]:
    """The CRAM of each sample that has one, as paired alignments."""
    tracks = []
    for sid in state.samples():
        url = _sample_file(state, sid, f"{sid}.GRCh38_GIABv3.cram")
        if url:
            tracks.append(
                {
                    "name": state.label(sid),
                    "type": "alignment",
                    "format": "cram",
                    "url": url,
                    "indexURL": f"{url}.crai",
                    "height": 250,
                    "displayMode": "SQUISHED",
                    "viewAsPairs": True,
                    "visibilityWindow": visibility_window_bp(),
                    "autoscale": False,
                    "autoscaleGroup": "cram",
                }
            )
    return tracks


def _resize_js(browser_var: str) -> str:
    return f"""
        setTimeout(() => {{
            if (window.{browser_var} && window.{browser_var}.visibilityChange) {{
                window.{browser_var}.visibilityChange();
            }}
        }}, 100);
    """


class IgvView:
    """One igv.js browser in a div, rebuilt by ``refresh``.

    The div is sized for the samples shown: 200 px, ``track_height`` per
    sample, 50 px of buffer.
    """

    def __init__(
        self,
        state: SVViewState,
        *,
        div_id: str,
        browser_var: str,
        tracks: Callable[[SVViewState], List[Dict[str, Any]]],
        track_height: int,
        locus: Callable[[], str],
        options: str = "",
    ) -> None:
        self.state = state
        self.div_id = div_id
        self.browser_var = browser_var
        self.tracks = tracks
        self.track_height = track_height
        self.locus = locus
        self.options = options
        self.container = (
            ui.element("div")
            .props(f'id="{div_id}"')
            .classes("w-full border border-gray-300 rounded-lg")
            .style(f"height: {self._height()}px")
        )

    def _height(self) -> int:
        return 200 + len(self.state.samples()) * self.track_height + 50

    def refresh(self) -> None:
        """Recreate the browser with the current samples' tracks."""
        self.container.clear()
        try:
            self.container.style(f"height: {self._height()}px")
        except RuntimeError:
            pass
        tracks = self.tracks(self.state)
        locus = self.locus()
        state = self.state
        roi = _roi_js(
            state, f"{state.start}-{state.end}", state.roi_start, state.roi_end
        )
        ui.run_javascript(
            f"""
            (async function() {{
                const igvDiv = document.getElementById("{self.div_id}");
                if (!igvDiv) return;

                // Clear any existing IGV instance
                if (window.{self.browser_var}) {{
                    igv.removeBrowser(window.{self.browser_var});
                }}

                const options = {{
                    genome: "hg38",
                    locus: "{locus}",
                    {self.options}
                    tracks: {json.dumps(tracks)},
                    roi: {roi},
                }};

                window.{self.browser_var} = await igv.createBrowser(igvDiv, options);

                // Trigger resize to ensure proper width
                {_resize_js(self.browser_var)}
            }})();
            """,
            timeout=10.0,
        )

    def resize(self) -> None:
        """Let igv.js measure its div again, after its panel was opened."""
        ui.run_javascript(_resize_js(self.browser_var))


def coverage_view(state: SVViewState) -> IgvView:
    """The by1000 bedgraph CNV view, the fallback without a coverage bigWig."""
    return IgvView(
        state,
        div_id="igv-div",
        browser_var="igvBrowser",
        tracks=coverage_tracks,
        track_height=100,
        locus=lambda: state.igv_locus,
    )


def vaf_view(state: SVViewState) -> IgvView:
    return IgvView(
        state,
        div_id="igv-vaf-div",
        browser_var="igvVafBrowser",
        tracks=vaf_tracks,
        track_height=100,
        locus=lambda: state.igv_locus,
    )


def _locus_center(locus: str) -> Optional[int]:
    """The middle of an igv.js ``chr:start-end`` locus, or None."""
    parts = locus.split(":")
    if len(parts) != 2:
        return None
    range_parts = parts[1].split("-")
    if len(range_parts) != 2:
        return None
    try:
        return round((float(range_parts[0]) + float(range_parts[1])) / 2)
    except ValueError:
        return None


class ReadLevelView:
    """The CRAM reads in the read windows, with New Start / New End pickers.

    A picker takes the centre of its window (the first for the start, the
    second for the end) and hands it to ``on_pick(key, position)``.
    """

    def __init__(self, state: SVViewState, on_pick: Callable[[str, int], None]) -> None:
        self.state = state
        self.on_pick = on_pick
        self.mode_label = ui.label().classes("text-sm text-gray-500 px-2 pb-2")
        # Curated coordinates when available; a space-separated locus list is
        # igv.js's multi-locus (split) view.
        self.view = IgvView(
            state,
            div_id="igv-cram-div",
            browser_var="igvCramBrowser",
            tracks=cram_tracks,
            track_height=250,
            locus=self._locus,
            options="showCenterGuide: true,",
        )
        self.refresh()
        with ui.row().classes("w-full items-center justify-center gap-0 mt-4"):
            ui.button(
                "New Start", on_click=lambda: self._pick("start", 0, "New Start")
            ).props("color=primary")
            ui.element("div").classes("w-[40%]")  # 40% spacer
            ui.button(
                "New End", on_click=lambda: self._pick("end", 1, "New End")
            ).props("color=primary")

    def _locus(self) -> str:
        used_mode, windows = self.state.read_windows()
        self.mode_label.text = self.state.describe_read_mode(used_mode)
        return " ".join(w.locus for w in windows)

    def refresh(self) -> None:
        self.view.refresh()

    def resize(self) -> None:
        self.view.resize()

    async def _pick(self, key: str, index: int, label: str) -> None:
        loci = await ui.run_javascript(
            "return window.igvCramBrowser ? window.igvCramBrowser.currentLoci() : null;",
            timeout=2.0,
        )
        if not (loci and isinstance(loci, list) and len(loci) > index):
            ui.notify("IGV browser not initialized", type="warning")
            return
        locus_str = loci[index]
        position = _locus_center(locus_str)
        if position is None:
            ui.notify(f"{label} - {locus_str}", type="info")
            return
        self.on_pick(key, position)
        ui.notify(f"{label} - {position}", type="info")


def update_roi_js(state: SVViewState) -> str:
    """Move the ROI of the CNV and read views to the current bounds.

    The read view is also re-centred on the read windows of those bounds.
    """
    bounds = state.bounds
    roi = _roi_js(
        state, f"{bounds['start']}-{bounds['end']}", bounds["start"], bounds["end"]
    )
    read_locus = " ".join(w.locus for w in state.read_windows()[1])
    return f"""
    const roiConfig = {roi};
    if (window.igvBrowser) {{
        window.igvBrowser.clearROIs();
        window.igvBrowser.loadROI(roiConfig);
    }}
    if (window.igvCramBrowser) {{
        window.igvCramBrowser.clearROIs();
        window.igvCramBrowser.loadROI(roiConfig);
        // Re-center the read windows on the new boundaries
        window.igvCramBrowser.search("{read_locus}");
    }}
    """
