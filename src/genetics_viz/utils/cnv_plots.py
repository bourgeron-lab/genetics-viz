"""ECharts option builder for the SV dialog's copy-number panel.

A pure function of the :class:`~genetics_viz.utils.cnv_coverage.CNTrack` list,
so the tests assert on the option dict rather than on a rendered chart.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

from genetics_viz.utils.cnv_coverage import CN_CONFIG, CNTrack, Panel

#: Pixel layout of one sample's grid.
GRID_HEIGHT = 130
GRID_GAP = 46
TOP = 34
BOTTOM = 64


def _mb_formatter(span_bp: int) -> str:
    """Tick labels in Mb, with enough decimals to tell the ticks apart."""
    decimals = 3 if span_bp < 1_000_000 else 2 if span_bp < 10_000_000 else 1
    return f"v => (v / 1e6).toFixed({decimals}) + ' Mb'"


def cn_chart_height(n_samples: int) -> int:
    return TOP + n_samples * (GRID_HEIGHT + GRID_GAP) + BOTTOM


def sample_subtitle(
    track: CNTrack, label: str, panel: Optional[Panel]
) -> Dict[str, Any]:
    """Title line of one grid: who, karyotype, fit, and why not divided."""
    parts = [f"{track.sample_id} {label}".strip()]
    if track.karyotype:
        parts.append(track.karyotype)
    rich_key = "ok"
    if track.fit_mad is not None:
        parts.append(f"fit MAD {track.fit_mad:.3f}")
        if (
            panel
            and panel.max_fit_mad is not None
            and track.fit_mad > panel.max_fit_mad
        ):
            rich_key = "warn"
    if track.mode == "corrected" and panel is not None:
        parts.append(f"÷ panel {panel.label}")
    text = " · ".join(parts)
    if track.reason:
        text = f"{text}{{warn|  ⚠ {track.reason}}}"
    return {"text": text, "rich_key": rich_key}


def cn_option(
    tracks: Sequence[CNTrack],
    labels: Sequence[str],
    *,
    window: tuple[int, int],
    roi: tuple[int, int],
    call: str,
    panel: Optional[Panel],
    cfg: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """One grid per sample, stacked, sharing the x-axis and its zoom."""
    cfg = cfg or CN_CONFIG
    colors = cfg.get("colors", {})
    y_max = float(cfg.get("y_max", 5))
    is_gain = "GAIN" in str(call).upper() or str(call).lower() == "dup"
    roi_color = colors.get("roi_gain" if is_gain else "roi_loss", "rgba(255,0,0,0.12)")
    n = len(tracks)
    mb_formatter = _mb_formatter(window[1] - window[0])

    grids: List[Dict[str, Any]] = []
    x_axes: List[Dict[str, Any]] = []
    y_axes: List[Dict[str, Any]] = []
    titles: List[Dict[str, Any]] = []
    series: List[Dict[str, Any]] = []

    for i, (track, label) in enumerate(zip(tracks, labels)):
        top = TOP + i * (GRID_HEIGHT + GRID_GAP)
        last = i == n - 1
        grids.append({"left": 56, "right": 24, "top": top, "height": GRID_HEIGHT})
        sub = sample_subtitle(track, label, panel)
        titles.append(
            {
                "text": sub["text"],
                "left": 56,
                "top": top - 22,
                "textStyle": {
                    "fontSize": 12,
                    "fontWeight": "normal",
                    "color": "#b45309" if sub["rich_key"] == "warn" else "#111827",
                    "rich": {"warn": {"color": "#b45309", "fontSize": 11}},
                },
            }
        )
        x_axes.append(
            {
                "type": "value",
                "gridIndex": i,
                "min": window[0],
                "max": window[1],
                "axisLabel": {"show": last, ":formatter": mb_formatter},
                "axisTick": {"show": last},
                "splitLine": {"show": False},
            }
        )
        y_axes.append(
            {
                "type": "value",
                "gridIndex": i,
                "min": 0,
                "max": y_max,
                "interval": 1,
                "name": "CN" if track.mode == "corrected" else "2×depth/med",
                "nameLocation": "middle",
                "nameGap": 30,
                "nameTextStyle": {"fontSize": 10, "color": "#6b7280"},
                "splitLine": {"lineStyle": {"color": "#f3f4f6"}},
            }
        )
        good = [
            [p, v]
            for p, v, m in zip(track.positions, track.values, track.masked)
            if v is not None and not m
        ]
        masked = [
            [p, v]
            for p, v, m in zip(track.positions, track.values, track.masked)
            if v is not None and m
        ]
        series.append(
            {
                "type": "scatter",
                "name": f"{track.sample_id}",
                "xAxisIndex": i,
                "yAxisIndex": i,
                "symbolSize": 3,
                "itemStyle": {"color": colors.get("neutral", "#4b5563")},
                "data": good,
                "markArea": {
                    "silent": True,
                    "itemStyle": {"color": roi_color},
                    "data": [[{"xAxis": roi[0]}, {"xAxis": roi[1]}]],
                },
                "markLine": {
                    "silent": True,
                    "symbol": "none",
                    "label": {"show": False},
                    "lineStyle": {"type": "dashed", "color": "#9ca3af", "width": 1},
                    "data": [{"yAxis": 1}, {"yAxis": 2}, {"yAxis": 3}],
                },
            }
        )
        series.append(
            {
                "type": "scatter",
                "name": f"{track.sample_id} (masked)",
                "xAxisIndex": i,
                "yAxisIndex": i,
                "symbolSize": 3,
                "itemStyle": {"color": colors.get("masked", "#d1d5db")},
                "data": masked,
            }
        )

    all_axes = list(range(n))
    return {
        "animation": False,
        "title": titles,
        "tooltip": {
            "trigger": "item",
            ":formatter": (
                "p => p.seriesName + '<br>' + (p.value[0] / 1e6).toFixed(4)"
                " + ' Mb: ' + p.value[1].toFixed(2)"
            ),
        },
        "grid": grids,
        "xAxis": x_axes,
        "yAxis": y_axes,
        "dataZoom": [
            {"type": "inside", "xAxisIndex": all_axes, "filterMode": "none"},
            {
                "type": "slider",
                "xAxisIndex": all_axes,
                "filterMode": "none",
                "bottom": 8,
                "height": 16,
                ":labelFormatter": mb_formatter,
            },
        ],
        "series": series,
    }
