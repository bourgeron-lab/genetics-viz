"""ECharts option builder for the SV dialog's read-evidence (samplot-like) panel.

The chart has one row per sample and one column per read window (one, or one
per breakpoint). Each cell is a single grid with two y-axes:

- **Left, insert size (log):**
  - each discordant pair as a thin line between its two read bodies;
  - each split read as its two aligned parts joined by an arc;
  - soft-clip clusters as bars at the foot;
  - interchromosomal mates as triangles at the top;
  - the expected band of concordant inserts, shaded.
- **Right, copy number:** read from the normalised coverage bigWig, divided by
  the cohort's control panel (``utils/cnv_coverage.py``).

A pair that crosses from one breakpoint window to the other runs off the inner
edge of each, as samplot's ``--zoom`` draws it.

Pure: tests assert on the returned dict.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, NamedTuple, Optional, Sequence, Tuple

from genetics_viz.utils.cnv_coverage import CN_CONFIG, CNTrack, Panel
from genetics_viz.utils.sv_reads import READ_CONFIG, WindowReads

SERIES_NAMES = {
    "deletion": "Deletion (+- large insert)",
    "duplication": "Duplication (-+)",
    "inversion_pp": "Inversion (++)",
    "inversion_mm": "Inversion (--)",
    "split": "Split read",
    "split_inv": "Split read (inverted)",
    "soft_clip": "Soft clips",
    "interchromosomal": "Interchromosomal mate",
    "insert_band": "Expected insert",
    "copy_number": "Copy number",
    "read_depth": "Read depth (2×depth/median)",
    "suggested": "Suggested breakpoint",
}

READ_HEIGHT = 210
ROW_GAP = 52
TOP = 84
BOTTOM = 58
LEFT_PCT = 6.0
RIGHT_PCT = 5.0
COL_GAP_PCT = 3.0
INSERT_FLOOR = 50
ARC_POINTS = 32
#: Evidence from reads under min_mapq, drawn when asked for.
LOW_MAPQ_OPACITY = 0.22

_POS_FORMATTER = "v => (v / 1e6).toFixed(4)"
_SIZE_FORMATTER = (
    "v => v >= 1e6 ? (v / 1e6) + ' Mb' : v >= 1e3 ? (v / 1e3) + ' kb' : v + ' bp'"
)
_TOOLTIP = (
    "p => (p.data && p.data.tip) ? p.data.tip"
    " : p.seriesName + ': ' + (Array.isArray(p.value) ? p.value[1] : p.value)"
)
#: Pan by dragging, no wheel zoom: the wheel scrolls the dialog.
INSIDE_ZOOM = {
    "type": "inside",
    "filterMode": "none",
    "zoomOnMouseWheel": False,
    "moveOnMouseWheel": False,
    "moveOnMouseMove": True,
}


def read_chart_height(n_samples: int) -> int:
    """No row gap after the last row: the axis labels and rulers follow it."""
    return TOP + n_samples * (READ_HEIGHT + ROW_GAP) - ROW_GAP + BOTTOM


def _fmt(n: float) -> str:
    return f"{int(round(n)):,}"


def _arc(x1: float, x2: float, y: float, rise: float) -> List[List[float]]:
    """Half an ellipse from ``x1`` to ``x2`` in log-y space, ``rise``x above ``y``.

    Steep at both ends, so it leaves the read parts almost vertically. The
    points crowd at the ends (cosine spacing), where the curve turns.
    """
    points = []
    for i in range(ARC_POINTS + 1):
        t = (1 - math.cos(math.pi * i / ARC_POINTS)) / 2
        height = math.sqrt(max(0.0, 1 - (2 * t - 1) ** 2))
        points.append([round(x1 + (x2 - x1) * t), round(y * rise**height, 1)])
    return points


class Suggested(NamedTuple):
    pos: int
    approximate: bool
    label: str


class Cell(NamedTuple):
    """One sample in one window."""

    reads: Optional[WindowReads]
    cn: Optional[CNTrack]
    #: ``[pos, 2 x depth / autosomal median]`` from the reads, unpanelled.
    depth: Optional[List[List[float]]] = None


def _row_title(label: str, cells: Sequence[Cell], panel: Optional[Panel]) -> str:
    parts = [label]
    reads = [c.reads for c in cells if c.reads is not None]
    cns = [c.cn for c in cells if c.cn is not None]
    if not reads:
        parts.append("no CRAM")
    else:
        med = next((r.insert_median for r in reads if r.insert_median), None)
        if med:
            parts.append(f"insert median {int(med)} bp")
    if cns:
        track = cns[0]
        if track.mode == "corrected" and panel is not None:
            parts.append(f"CN ÷ panel {panel.label}")
        elif track.reason:
            parts.append(f"⚠ {track.reason}")
    else:
        parts.append("no coverage file")
    return "  ·  ".join(parts)


def read_option(
    rows: Sequence[Sequence[Cell]],
    labels: Sequence[str],
    *,
    windows: Sequence[Any],
    breakpoints: Sequence[int],
    suggested: Sequence[Suggested] = (),
    intervals: Sequence[Tuple[int, int]] = (),
    panel: Optional[Panel] = None,
    cfg: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Grids, axes and series for every sample x window cell.

    ``breakpoints`` are the SV's current bounds (curated when set), drawn
    dotted red. ``suggested`` are the suggested breakpoints, drawn green:
    solid when exact (reads), dashed when approximate (coverage step).
    ``intervals`` are the ranges the discordant pairs allow, shaded green.
    """
    rows = [[Cell(*c) for c in cells] for cells in rows]
    cfg = cfg or READ_CONFIG
    colors = cfg.get("colors", {})
    rise = float(cfg.get("arc_rise", 1.8))
    cn_max = float(CN_CONFIG.get("y_max", 5))
    n_windows = len(windows)
    col_width = (100 - LEFT_PCT - RIGHT_PCT - COL_GAP_PCT * (n_windows - 1)) / max(
        1, n_windows
    )

    def col_left(col: int) -> float:
        return LEFT_PCT + col * (col_width + COL_GAP_PCT)

    grids: List[Dict[str, Any]] = []
    x_axes: List[Dict[str, Any]] = []
    y_axes: List[Dict[str, Any]] = []
    series: List[Dict[str, Any]] = []
    titles: List[Dict[str, Any]] = [
        {
            "text": w.locus,
            "left": f"{col_left(col) + col_width / 2}%",
            "textAlign": "center",
            "top": 6,
            "textStyle": {"fontSize": 12, "color": "#374151"},
        }
        for col, w in enumerate(windows)
    ]
    columns: List[List[int]] = [[] for _ in range(n_windows)]

    for row, (cells, label) in enumerate(zip(rows, labels)):
        top = TOP + row * (READ_HEIGHT + ROW_GAP)
        heights = [
            y
            for c in cells
            if c.reads is not None
            for y in [p.y for p in c.reads.pairs] + [s.y * rise for s in c.reads.splits]
        ]
        y_top = max(heights + [2000]) * 1.5
        titles.append(
            {
                "text": _row_title(label, cells, panel),
                "left": f"{LEFT_PCT}%",
                "top": top - 22,
                "textStyle": {"fontSize": 12, "fontWeight": "normal"},
            }
        )
        last_row = row == len(rows) - 1
        for col, (cell, w) in enumerate(zip(cells, windows)):
            idx = len(grids)
            ins_axis, cn_axis = len(y_axes), len(y_axes) + 1
            last_col = col == n_windows - 1
            grids.append(
                {
                    "left": f"{col_left(col)}%",
                    "width": f"{col_width}%",
                    "top": top,
                    "height": READ_HEIGHT,
                }
            )
            x_axes.append(
                {
                    "type": "value",
                    "gridIndex": idx,
                    "min": w.start,
                    "max": w.end,
                    # The end labels of side-by-side windows would collide.
                    "axisLabel": {
                        "show": last_row,
                        "fontSize": 10,
                        "showMinLabel": False,
                        "showMaxLabel": False,
                        ":formatter": _POS_FORMATTER,
                    },
                    "axisTick": {"show": last_row},
                    "splitLine": {"show": False},
                }
            )
            columns[col].append(idx)
            y_axes.append(
                {
                    "type": "log",
                    "gridIndex": idx,
                    "min": INSERT_FLOOR,
                    # A power of ten, so the top label is a round size.
                    "max": 10 ** math.ceil(math.log10(y_top * 1.3)),
                    "name": "insert size" if col == 0 else "",
                    "nameLocation": "middle",
                    "nameGap": 40,
                    "nameTextStyle": {"fontSize": 10, "color": "#6b7280"},
                    "axisLabel": {
                        "show": col == 0,
                        "fontSize": 9,
                        ":formatter": _SIZE_FORMATTER,
                    },
                    "splitLine": {"lineStyle": {"color": "#f3f4f6"}},
                }
            )
            y_axes.append(
                {
                    "type": "value",
                    "gridIndex": idx,
                    "position": "right",
                    "min": 0,
                    "max": cn_max,
                    "interval": 1,
                    "name": "CN" if last_col else "",
                    "nameLocation": "middle",
                    "nameGap": 22,
                    "nameTextStyle": {
                        "fontSize": 10,
                        "color": colors.get("copy_number"),
                    },
                    "axisLabel": {
                        "show": last_col,
                        "fontSize": 9,
                        "color": colors.get("copy_number"),
                    },
                    "splitLine": {"show": False},
                }
            )
            series += _cell_series(
                cell,
                w,
                breakpoints,
                suggested,
                intervals,
                ins_axis=ins_axis,
                cn_axis=cn_axis,
                x_axis=idx,
                y_top=y_top,
                rise=rise,
                colors=colors,
            )

    present: List[str] = []
    for s in series:
        if s["name"] not in present and (s.get("data") or s.get("markArea")):
            present.append(s["name"])
    legend = [n for n in SERIES_NAMES.values() if n in present]
    return {
        "animation": False,
        "title": titles,
        # One line whatever the width, so it never runs into the first row.
        "legend": {
            "type": "scroll",
            "data": legend,
            "top": 30,
            "left": f"{LEFT_PCT}%",
            "right": f"{RIGHT_PCT}%",
            "itemWidth": 18,
            "textStyle": {"fontSize": 11},
        },
        "tooltip": {"trigger": "item", ":formatter": _TOOLTIP},
        "grid": grids,
        "xAxis": x_axes,
        "yAxis": y_axes,
        "dataZoom": [{**INSIDE_ZOOM, "xAxisIndex": axes} for axes in columns]
        + [
            {
                "type": "slider",
                "xAxisIndex": axes,
                "filterMode": "none",
                "bottom": 10,
                "height": 16,
                "left": f"{col_left(col)}%",
                "width": f"{col_width}%",
                ":labelFormatter": _POS_FORMATTER,
            }
            for col, axes in enumerate(columns)
        ],
        "series": series,
    }


def _cell_series(
    cell: Cell,
    window,
    breakpoints: Sequence[int],
    suggested: Sequence[Suggested],
    intervals: Sequence[Tuple[int, int]],
    *,
    ins_axis: int,
    cn_axis: int,
    x_axis: int,
    y_top: float,
    rise: float,
    colors: Dict[str, str],
) -> List[Dict[str, Any]]:
    reads, cn = cell.reads, cell.cn
    out: List[Dict[str, Any]] = []
    on_ins = {"xAxisIndex": x_axis, "yAxisIndex": ins_axis}
    on_cn = {"xAxisIndex": x_axis, "yAxisIndex": cn_axis}

    # The band and the breakpoints ride on an empty series of the insert axis.
    band: Dict[str, Any] = {
        "type": "line",
        "name": SERIES_NAMES["insert_band"],
        **on_ins,
        "data": [],
        "markLine": {
            "silent": True,
            "symbol": "none",
            "label": {"show": False},
            "lineStyle": {"color": "#ef4444", "type": "dotted", "width": 1},
            "data": [
                {"xAxis": b} for b in breakpoints if window.start <= b <= window.end
            ]
            + [
                {
                    "xAxis": sug.pos,
                    "lineStyle": {
                        "color": colors.get("suggested_breakpoint", "#16a34a"),
                        "type": "dashed" if sug.approximate else "solid",
                        "width": 1.5,
                    },
                    "label": {
                        "show": True,
                        "position": "insideEndTop",
                        "formatter": sug.label,
                        "fontSize": 10,
                        "color": colors.get("suggested_breakpoint", "#16a34a"),
                    },
                }
                for sug in suggested
                if window.start <= sug.pos <= window.end
            ],
        },
    }
    areas: List[Any] = []
    if reads is not None and reads.insert_lo and reads.insert_hi:
        areas.append([{"yAxis": reads.insert_lo}, {"yAxis": reads.insert_hi}])
    for lo, hi in intervals:
        if hi >= window.start and lo <= window.end:
            areas.append(
                [
                    {
                        "xAxis": max(lo, window.start),
                        "itemStyle": {
                            "color": colors.get(
                                "breakpoint_interval", "rgba(22,163,74,0.1)"
                            )
                        },
                    },
                    {"xAxis": min(hi, window.end)},
                ]
            )
    if areas:
        band["markArea"] = {
            "silent": True,
            "itemStyle": {"color": colors.get("insert_band", "rgba(0,0,0,0.08)")},
            "data": areas,
        }
    out.append(band)

    if cn is not None:
        line = {
            "type": "line",
            "step": "middle",
            "showSymbol": False,
            "connectNulls": False,
            "z": 1,
            **on_cn,
        }
        out.append(
            {
                **line,
                "name": SERIES_NAMES["copy_number"],
                "lineStyle": {"width": 1.5, "color": colors.get("copy_number")},
                "itemStyle": {"color": colors.get("copy_number")},
                "areaStyle": {"color": colors.get("copy_number"), "opacity": 0.08},
                "data": [
                    [p, None if m else v]
                    for p, v, m in zip(cn.positions, cn.values, cn.masked)
                ],
                "markLine": {
                    "silent": True,
                    "symbol": "none",
                    "label": {"show": False},
                    "lineStyle": {"color": "#93c5fd", "type": "dashed", "width": 1},
                    "data": [{"yAxis": 2}],
                },
            }
        )
        if any(cn.masked):
            out.append(
                {
                    **line,
                    "name": SERIES_NAMES["copy_number"],
                    "lineStyle": {
                        "width": 1.5,
                        "color": colors.get("copy_number_masked"),
                    },
                    "data": [
                        [p, v if m else None]
                        for p, v, m in zip(cn.positions, cn.values, cn.masked)
                    ],
                }
            )

    if cell.depth:
        cn_max = float(CN_CONFIG.get("y_max", 5))
        out.append(
            {
                "type": "line",
                "name": SERIES_NAMES["read_depth"],
                **on_cn,
                "showSymbol": False,
                "z": 2,
                "lineStyle": {
                    "width": 1,
                    "color": colors.get("read_depth", "#0f766e"),
                    "opacity": 0.8,
                },
                "itemStyle": {"color": colors.get("read_depth", "#0f766e")},
                "data": [[p, min(v, cn_max)] for p, v in cell.depth],
            }
        )

    if reads is None:
        return out

    lines = {"type": "lines", "coordinateSystem": "cartesian2d", **on_ins, "z": 3}
    for cls, low in (
        (c, low)
        for c in ("deletion", "duplication", "inversion_pp", "inversion_mm")
        for low in (False, True)
    ):
        pairs = [p for p in reads.pairs if p.cls == cls and p.low == low]
        if not pairs:
            continue
        color = colors.get(cls, "#000")
        name = SERIES_NAMES[cls]
        fade = LOW_MAPQ_OPACITY if low else 1.0
        out.append(
            {
                **lines,
                "name": name,
                "lineStyle": {"color": color, "width": 1, "opacity": 0.8 * fade},
                "data": [
                    {
                        "coords": [[p.x1, p.y], [p.x2, p.y]],
                        "tip": f"{name}<br>insert {_fmt(p.y)} bp",
                    }
                    for p in pairs
                ],
            }
        )
        out.append(
            {
                **lines,
                "name": name,
                "lineStyle": {"color": color, "width": 3.5, "opacity": 0.95 * fade},
                "data": [
                    {
                        "coords": [[a, p.y], [b, p.y]],
                        "tip": f"{name}<br>insert {_fmt(p.y)} bp",
                    }
                    for p in pairs
                    for a, b in ((p.x1, p.left_end), (p.right_start, p.x2))
                ],
            }
        )

    for cls, low in ((c, low) for c in ("split", "split_inv") for low in (False, True)):
        splits = [s for s in reads.splits if s.cls == cls and s.low == low]
        if not splits:
            continue
        color = colors.get(cls, "#000")
        name = SERIES_NAMES[cls]
        fade = LOW_MAPQ_OPACITY if low else 1.0
        out.append(
            {
                **lines,
                "name": name,
                "z": 4,
                "lineStyle": {"color": color, "width": 3.5, "opacity": 0.9 * fade},
                "data": [
                    {
                        "coords": [[a, s.y], [b, s.y]],
                        "tip": f"{name}<br>gap {_fmt(s.y)} bp",
                    }
                    for s in splits
                    for a, b in ((s.a1, s.a2), (s.b1, s.b2))
                ],
            }
        )
        # The link arcs above the parts, so it never hides a pair's line.
        out.append(
            {
                **lines,
                "name": name,
                "z": 2,
                "polyline": True,
                "lineStyle": {"color": color, "width": 1, "opacity": 0.7 * fade},
                "data": [
                    {
                        "coords": _arc(s.a2, s.b1, s.y, rise),
                        "tip": f"{name}<br>gap {_fmt(s.y)} bp",
                    }
                    for s in splits
                ],
            }
        )

    if reads.clips:
        biggest = max(c.count for c in reads.clips)
        out.append(
            {
                "type": "scatter",
                "name": SERIES_NAMES["soft_clip"],
                **on_ins,
                "symbol": "rect",
                "z": 5,
                "itemStyle": {"color": colors.get("soft_clip", "#f59e0b")},
                "data": [
                    {
                        "value": [c.pos, INSERT_FLOOR * 1.25],
                        "itemStyle": {"opacity": LOW_MAPQ_OPACITY if c.low else 1.0},
                        "symbolSize": [3, round(6 + 18 * c.count / biggest)],
                        "tip": (
                            f"{c.count} reads soft-clipped at {_fmt(c.pos)}"
                            f" ({'before' if c.side == 'left' else 'after'} it)"
                            + (", MAPQ < min" if c.low else "")
                        ),
                    }
                    for c in reads.clips
                ],
            }
        )

    if reads.interchrom:
        out.append(
            {
                "type": "scatter",
                "name": SERIES_NAMES["interchromosomal"],
                **on_ins,
                "symbol": "triangle",
                "symbolSize": 7,
                "itemStyle": {"color": colors.get("interchromosomal")},
                "data": [
                    {"value": [m.x, y_top], "tip": f"Mate on {m.mate}"}
                    for m in reads.interchrom
                ],
            }
        )
    return out
