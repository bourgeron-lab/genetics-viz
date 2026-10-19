"""ECharts option builder for the SV dialog's read-evidence (samplot-like) panel.

One row per sample, one column per read window (one, or one per breakpoint).
Each cell is a read grid (pairs and split reads as horizontal segments, insert
size on a log y-axis) over a small coverage grid. A pair crossing from one
breakpoint window to the other runs off the inner edge of each, which is what
samplot's ``--zoom`` draws too.

Pure: tests assert on the returned dict.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

from genetics_viz.utils.sv_reads import READ_CONFIG, WindowReads

SERIES_NAMES = {
    "concordant": "Concordant",
    "deletion": "Deletion (+- large insert)",
    "duplication": "Duplication (-+)",
    "inversion_pp": "Inversion (++)",
    "inversion_mm": "Inversion (--)",
    "split": "Split read",
    "split_inv": "Split read (inverted)",
    "interchromosomal": "Interchromosomal mate",
    "coverage": "Coverage",
}
_COLOR_KEYS = {"split_inv": "inversion_pp"}

READ_HEIGHT = 150
COV_HEIGHT = 36
ROW_GAP = 44
TOP = 80
BOTTOM = 20
LEFT_PCT = 6.0
RIGHT_PCT = 2.0
COL_GAP_PCT = 4.0

_KB_FORMATTER = "v => (v / 1e6).toFixed(4)"
_SIZE_FORMATTER = (
    "v => v >= 1e6 ? (v / 1e6) + ' Mb' : v >= 1e3 ? (v / 1e3) + ' kb' : v + ' bp'"
)
_TOOLTIP = (
    "p => { if (p.seriesType === 'lines') return p.seriesName + '<br>'"
    " + p.data.coords[0][0] + '-' + p.data.coords[1][0] + ', ' + p.data.coords[0][1]"
    " + ' bp'; if (p.data && p.data.mate) return 'Mate on ' + p.data.mate;"
    " return p.seriesName + ': ' + p.value[1]; }"
)


def read_chart_height(n_samples: int) -> int:
    return TOP + n_samples * (READ_HEIGHT + COV_HEIGHT + ROW_GAP) + BOTTOM


def _color(cls: str, colors: Dict[str, str]) -> str:
    return colors.get(cls) or colors.get(_COLOR_KEYS.get(cls, ""), "#000000")


def read_option(
    samples: Sequence[List[WindowReads]],
    labels: Sequence[str],
    *,
    breakpoints: Sequence[int],
    cfg: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Grids, axes and series for every sample x window cell."""
    cfg = cfg or READ_CONFIG
    colors = cfg.get("colors", {})
    n_windows = len(samples[0]) if samples else 0
    col_width = (100 - LEFT_PCT - RIGHT_PCT - COL_GAP_PCT * (n_windows - 1)) / max(
        1, n_windows
    )

    grids: List[Dict[str, Any]] = []
    x_axes: List[Dict[str, Any]] = []
    y_axes: List[Dict[str, Any]] = []
    series: List[Dict[str, Any]] = []
    titles: List[Dict[str, Any]] = []
    zoom_axes: List[List[int]] = [[] for _ in range(n_windows)]

    for row, (cells, label) in enumerate(zip(samples, labels)):
        top = TOP + row * (READ_HEIGHT + COV_HEIGHT + ROW_GAP)
        highest = max(
            [s.y for c in cells for s in (*c.pairs, *c.splits)] + [1000],
        )
        y_top = highest * 1.6
        ins = [c.insert_median for c in cells if c.insert_median]
        sub = f"insert median {int(ins[0])} bp" if ins else ""
        titles.append(
            {
                "text": f"{label}" + (f"  ·  {sub}" if sub else ""),
                "left": f"{LEFT_PCT}%",
                "top": top - 22,
                "textStyle": {"fontSize": 12, "fontWeight": "normal"},
            }
        )
        last_row = row == len(samples) - 1
        for col, cell in enumerate(cells):
            left = LEFT_PCT + col * (col_width + COL_GAP_PCT)
            w = cell.window
            read_idx = len(grids)
            cov_idx = read_idx + 1
            grids.append(
                {
                    "left": f"{left}%",
                    "width": f"{col_width}%",
                    "top": top,
                    "height": READ_HEIGHT,
                }
            )
            grids.append(
                {
                    "left": f"{left}%",
                    "width": f"{col_width}%",
                    "top": top + READ_HEIGHT + 4,
                    "height": COV_HEIGHT,
                }
            )
            for idx, show in ((read_idx, False), (cov_idx, last_row)):
                x_axes.append(
                    {
                        "type": "value",
                        "gridIndex": idx,
                        "min": w.start,
                        "max": w.end,
                        # The end labels of side-by-side windows would collide.
                        "axisLabel": {
                            "show": show,
                            "fontSize": 10,
                            "showMinLabel": False,
                            "showMaxLabel": False,
                            ":formatter": _KB_FORMATTER,
                        },
                        "axisTick": {"show": show},
                        "splitLine": {"show": False},
                    }
                )
                zoom_axes[col].append(idx)
            y_axes.append(
                {
                    "type": "log",
                    "gridIndex": read_idx,
                    "min": 50,
                    "max": y_top * 1.3,
                    "name": "insert (bp)" if col == 0 else "",
                    "nameLocation": "middle",
                    "nameGap": 34,
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
                    "gridIndex": cov_idx,
                    "min": 0,
                    "axisLabel": {"show": False},
                    "splitLine": {"show": False},
                }
            )
            if row == 0:
                titles.append(
                    {
                        "text": w.locus,
                        "left": f"{left + col_width / 2}%",
                        "textAlign": "center",
                        "top": 6,
                        "textStyle": {"fontSize": 12, "color": "#374151"},
                    }
                )

            by_cls: Dict[str, List[Dict[str, Any]]] = {}
            for seg in (*cell.pairs, *cell.splits):
                by_cls.setdefault(seg.cls, []).append(
                    {"coords": [[seg.x1, seg.y], [seg.x2, seg.y]]}
                )
            for cls, data in by_cls.items():
                series.append(
                    {
                        "type": "lines",
                        "name": SERIES_NAMES[cls],
                        "coordinateSystem": "cartesian2d",
                        "xAxisIndex": read_idx,
                        "yAxisIndex": read_idx,
                        "polyline": False,
                        "lineStyle": {
                            "color": _color(cls, colors),
                            "width": 1 if cls == "concordant" else 1.5,
                            "opacity": 0.5 if cls == "concordant" else 0.9,
                            "type": "dashed" if cls.startswith("split") else "solid",
                        },
                        "z": 1 if cls == "concordant" else 3,
                        "data": data,
                    }
                )
            in_window = [b for b in breakpoints if w.start <= b <= w.end]
            series.append(
                {
                    "type": "scatter",
                    "name": SERIES_NAMES["interchromosomal"],
                    "xAxisIndex": read_idx,
                    "yAxisIndex": read_idx,
                    "symbol": "triangle",
                    "symbolSize": 7,
                    "itemStyle": {"color": _color("interchromosomal", colors)},
                    "data": [
                        {"value": [m.x, y_top], "mate": m.mate} for m in cell.interchrom
                    ],
                    "markLine": {
                        "silent": True,
                        "symbol": "none",
                        "label": {"show": False},
                        "lineStyle": {"color": "#ef4444", "type": "dotted", "width": 1},
                        "data": [{"xAxis": b} for b in in_window],
                    },
                }
            )
            series.append(
                {
                    "type": "line",
                    "name": SERIES_NAMES["coverage"],
                    "xAxisIndex": cov_idx,
                    "yAxisIndex": cov_idx,
                    "showSymbol": False,
                    "lineStyle": {"width": 0.5, "color": _color("coverage", colors)},
                    "areaStyle": {"color": _color("coverage", colors), "opacity": 0.5},
                    "data": cell.coverage,
                }
            )

    present = []
    for s in series:
        if s["name"] not in present:
            present.append(s["name"])
    legend = [n for n in SERIES_NAMES.values() if n in present]
    return {
        "animation": False,
        "title": titles,
        "legend": {
            "data": legend,
            "top": 28,
            "itemWidth": 18,
            "textStyle": {"fontSize": 11},
        },
        "tooltip": {"trigger": "item", ":formatter": _TOOLTIP},
        "grid": grids,
        "xAxis": x_axes,
        "yAxis": y_axes,
        "dataZoom": [
            {"type": "inside", "xAxisIndex": axes, "filterMode": "none"}
            for axes in zoom_axes
        ],
        "series": series,
    }
