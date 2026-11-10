"""Plots of the upd step: the chromosome view, the genome overview, the cohort chart.

Pure builders: they take table rows and return an ECharts option or an SVG
string. The counts come from the step's ``windows.tsv`` (``upd_v4`` and later):

- up, the alleles only the mother has (``n_mat``: father homozygous, child
  het) and, darker, the Mendelian errors against the father
  (``n_lacks_father``); down, the same for the father;
- the child's het rate, and its depth ratio on the right axis;
- lanes: the step's segments, the ROH step's runs, the Sample QC span, the
  validated span, and the cytobands.

A maternal UPD has no paternal-specific allele (iso or het); an isodisomy has
no het call and many errors against the other parent; a deletion of one copy
looks like an isodisomy, at half the depth.
"""

import html
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .cytobands import CYTOBANDS, GIESTAIN_COLORS
from .ideogram import chrom_length_mb, chromosome_svg, mb_axis_svg
from .upd_calls import (
    CONFIDENCES,
    WINDOW_STYLE,
    composition_label,
    confidence_info,
    origin_info,
    span_text,
)

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]

_MB = "v => (v / 1e6).toFixed(0) + ' Mb'"

#: Lane order, top to bottom (the y axis is a category axis, bottom first).
_LANES = ["Validated", "UPD segments", "ROH runs", "Sample QC span", "Bands"]


def _frac(value: Any, total: Any) -> Optional[float]:
    try:
        if value is None or not total:
            return None
        return round(float(value) / float(total), 4)
    except (TypeError, ValueError):
        return None


def _het_rate(w: Dict[str, Any]) -> Optional[float]:
    het, hom = w.get("n_het"), w.get("n_hom_alt")
    if het is None or hom is None or het + hom == 0:
        return None
    return round(het / (het + hom), 4)


def _with_gaps(windows: Sequence[Dict[str, Any]], value) -> List[Optional[List[float]]]:
    """``[[mid, y], …]`` with a null where a window is missing (masked genome)."""
    out: List[Optional[List[float]]] = []
    previous_end: Optional[int] = None
    for w in windows:
        start, end = int(w["start"]), int(w["end"])
        if previous_end is not None and start > previous_end + 1:
            out.append(None)
        y = value(w)
        out.append([(start + end) / 2, y] if y is not None else None)
        previous_end = end
    return out


def _lane(
    name: str, intervals: Sequence[Tuple[int, int]], color: str, label: str
) -> Dict[str, Any]:
    """One lane of thick horizontal strokes, ``null`` between intervals."""
    data: List[Optional[List[Any]]] = []
    for start, end in intervals:
        data += [[start, name], [end, name], None]
    return {
        "type": "line",
        "name": label,
        "xAxisIndex": 2,
        "yAxisIndex": 3,
        "symbol": "none",
        "connectNulls": False,
        "lineStyle": {"width": 10, "color": color, "cap": "butt"},
        "itemStyle": {"color": color},
        "data": data,
        "silent": True,
    }


def chromosome_option(
    chrom: str,
    windows: Sequence[Dict[str, Any]],
    *,
    segments: Sequence[Dict[str, Any]] = (),
    roh_runs: Sequence[Dict[str, Any]] = (),
    qc_span: Optional[Tuple[int, int]] = None,
    validated: Optional[Tuple[int, int]] = None,
    events: Sequence[Dict[str, Any]] = (),
) -> Dict[str, Any]:
    """The review dialog's view of one member's chromosome.

    Without ``windows`` (labels before ``upd_v4``) only the lanes are drawn.
    """
    style = WINDOW_STYLE
    length = chrom_length_mb(chrom) * 1e6
    has_windows = bool(windows)
    # Keep in step with chromosome_height.
    bars_top, bars_h = 40, 150
    lines_top, lines_h = bars_top + bars_h + 30, 90
    lanes_top = (lines_top + lines_h + 30) if has_windows else 40
    lanes_h = 26 * len(_LANES)

    grids = [
        {"left": 64, "right": 56, "top": bars_top, "height": bars_h},
        {"left": 64, "right": 56, "top": lines_top, "height": lines_h},
        {"left": 64, "right": 56, "top": lanes_top, "height": lanes_h},
    ]
    if not has_windows:
        grids[0] = {**grids[0], "show": False, "height": 1}
        grids[1] = {**grids[1], "show": False, "height": 1, "top": 20}
    x_axes = [
        {
            "type": "value",
            "gridIndex": i,
            "min": 0,
            "max": length,
            "show": i == 2 or has_windows,
            "axisLabel": {"show": i == 2, ":formatter": _MB},
            "axisTick": {"show": i == 2},
            "splitLine": {"show": False},
        }
        for i in range(3)
    ]
    y_axes = [
        {
            "type": "value",
            "gridIndex": 0,
            "show": has_windows,
            "name": "share of sites",
            "nameLocation": "middle",
            "nameGap": 44,
            "nameTextStyle": {"fontSize": 10, "color": "#6b7280"},
            "axisLabel": {":formatter": "v => Math.abs(v).toFixed(2)"},
            "splitLine": {"lineStyle": {"color": "#f3f4f6"}},
        },
        {
            "type": "value",
            "gridIndex": 1,
            "show": has_windows,
            "min": 0,
            "max": 1,
            "name": "het rate",
            "nameLocation": "middle",
            "nameGap": 44,
            "nameTextStyle": {"fontSize": 10, "color": "#6b7280"},
            "splitLine": {"lineStyle": {"color": "#f3f4f6"}},
        },
        {
            "type": "value",
            "gridIndex": 1,
            "show": has_windows,
            "position": "right",
            "min": 0,
            "max": 2,
            "name": "depth",
            "nameLocation": "middle",
            "nameGap": 36,
            "nameTextStyle": {"fontSize": 10, "color": style.get("depth_ratio")},
            "splitLine": {"show": False},
        },
        {
            "type": "category",
            "gridIndex": 2,
            "data": list(reversed(_LANES)),
            "axisTick": {"show": False},
            "axisLine": {"show": False},
            "axisLabel": {"fontSize": 10, "color": "#374151"},
        },
    ]

    series: List[Dict[str, Any]] = []
    if has_windows:
        shading = [
            [
                {"xAxis": e["start"], "itemStyle": {"color": _shade(e)}},
                {"xAxis": e["end"]},
            ]
            for e in events
            if e.get("start") is not None and e.get("end") is not None
        ]
        for name, column, color, sign in (
            ("Allele only from the mother", "n_mat", style.get("maternal"), 1),
            ("Lacks the father's allele", "n_lacks_father",
             style.get("maternal_error"), 1),
            ("Allele only from the father", "n_pat", style.get("paternal"), -1),
            ("Lacks the mother's allele", "n_lacks_mother",
             style.get("paternal_error"), -1),
        ):  # fmt: skip
            data = _with_gaps(
                windows,
                lambda w, c=column, s=sign: (
                    None if (f := _frac(w.get(c), w.get("n_sites"))) is None else s * f
                ),
            )
            if not any(p is not None for p in data):
                continue
            series.append(
                {
                    "type": "line",
                    "name": name,
                    "xAxisIndex": 0,
                    "yAxisIndex": 0,
                    "stack": "up" if sign > 0 else "down",
                    "step": "middle",
                    "symbol": "none",
                    "connectNulls": False,
                    "lineStyle": {"width": 0},
                    "areaStyle": {"color": color, "opacity": 0.9},
                    "itemStyle": {"color": color},
                    "data": data,
                }
            )
        if series:
            series[0]["markArea"] = {"silent": True, "data": shading}
            series[0]["markLine"] = {
                "silent": True,
                "symbol": "none",
                "label": {"show": False},
                "lineStyle": {"color": "#9ca3af", "width": 1},
                "data": [{"yAxis": 0}],
            }
        series.append(
            {
                "type": "line",
                "name": "Het rate",
                "xAxisIndex": 1,
                "yAxisIndex": 1,
                "step": "middle",
                "symbol": "none",
                "connectNulls": False,
                "lineStyle": {"width": 1.2, "color": style.get("het_rate")},
                "itemStyle": {"color": style.get("het_rate")},
                "data": _with_gaps(windows, _het_rate),
                "markArea": {"silent": True, "data": shading},
            }
        )
        series.append(
            {
                "type": "line",
                "name": "Depth ratio",
                "xAxisIndex": 1,
                "yAxisIndex": 2,
                "step": "middle",
                "symbol": "none",
                "connectNulls": False,
                "lineStyle": {"width": 1, "color": style.get("depth_ratio")},
                "itemStyle": {"color": style.get("depth_ratio")},
                "data": _with_gaps(windows, lambda w: w.get("depth_ratio")),
                "markLine": {
                    "silent": True,
                    "symbol": "none",
                    "label": {"show": False},
                    "lineStyle": {
                        "type": "dashed",
                        "color": style.get("depth_ratio"),
                        "width": 1,
                    },
                    "data": [
                        {"yAxis": 1},
                        {"yAxis": float(style.get("deletion_depth", 0.75))},
                    ],
                },
            }
        )

    # Lanes: one series per colour.
    by_kind: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for s in segments:
        if s.get("start") is None or s.get("end") is None:
            continue
        label = (
            f"{composition_label(s.get('state'))} "
            f"{origin_info(s.get('origin'))['label']}"
        )
        by_kind.setdefault((label, origin_info(s.get("origin"))["color"]), []).append(
            (int(s["start"]), int(s["end"]))
        )
    for (label, color), intervals in by_kind.items():
        series.append(_lane("UPD segments", intervals, color, f"Segments: {label}"))
    runs = [
        (int(r["start"]), int(r["end"]))
        for r in roh_runs
        if r.get("start") is not None and r.get("end") is not None
    ]
    if runs:
        series.append(_lane("ROH runs", runs, style.get("roh", "#6a1b9a"), "ROH runs"))
    if qc_span:
        series.append(
            _lane(
                "Sample QC span",
                [qc_span],
                style.get("qc_span", "#f9a825"),
                "Sample QC opposite homozygotes",
            )
        )
    if validated:
        series.append(_lane("Validated", [validated], "#2e7d32", "Validated span"))
    bands: Dict[str, List[Tuple[int, int]]] = {}
    for band in CYTOBANDS.get(chrom, []):
        bands.setdefault(band["stain"], []).append(
            (int(band["start"] * 1e6), int(band["end"] * 1e6))
        )
    for stain, intervals in bands.items():
        lane = _lane("Bands", intervals, GIESTAIN_COLORS.get(stain, "#e5e7eb"), stain)
        lane["lineStyle"]["width"] = 12
        lane["legendHoverLink"] = False
        series.append(lane)

    legend = [s["name"] for s in series if s["name"] not in bands]
    return {
        "animation": False,
        "legend": {
            "top": 0,
            "type": "scroll",
            "data": legend,
            "textStyle": {"fontSize": 11},
        },
        "tooltip": {
            "trigger": "axis",
            "axisPointer": {"type": "line"},
            ":formatter": (
                "ps => { const p = ps.filter(q => q.value && q.value[1] != null"
                " && typeof q.value[1] === 'number'); if (!p.length) return '';"
                " return (p[0].value[0] / 1e6).toFixed(2) + ' Mb<br>' + p.map(q =>"
                " q.marker + q.seriesName + ': ' + Math.abs(q.value[1]).toFixed(3))"
                ".join('<br>'); }"
            ),
        },
        "grid": grids,
        "xAxis": x_axes,
        "yAxis": y_axes,
        "dataZoom": [
            # Pan by dragging, no wheel zoom: the wheel scrolls the dialog.
            {
                "type": "inside",
                "xAxisIndex": [0, 1, 2],
                "filterMode": "none",
                "zoomOnMouseWheel": False,
                "moveOnMouseWheel": False,
                "moveOnMouseMove": True,
            },
            {
                "type": "slider",
                "xAxisIndex": [0, 1, 2],
                "filterMode": "none",
                "bottom": 8,
                "height": 16,
                ":labelFormatter": _MB,
            },
        ],
        "series": series,
    }


def chromosome_height(has_windows: bool) -> int:
    """Pixel height of :func:`chromosome_option`'s chart."""
    lanes = 26 * len(_LANES)
    return (40 + 150 + 30 + 90 + 30 if has_windows else 40) + lanes + 70


def _shade(event: Dict[str, Any]) -> str:
    color = confidence_info(event.get("confidence"))["color"]
    return _rgba(color, 0.08)


def _rgba(hex_color: str, alpha: float) -> str:
    h = hex_color.lstrip("#")
    if len(h) != 6:
        return f"rgba(158,158,158,{alpha})"
    r, g, b = (int(h[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


# ---------------------------------------------------------------------------
# Genome overview
# ---------------------------------------------------------------------------


def _mix(a: str, b: str, t: float) -> str:
    """``a`` at 0, ``b`` at 1."""
    a, b = a.lstrip("#"), b.lstrip("#")
    ca = [int(a[i : i + 2], 16) for i in (0, 2, 4)]
    cb = [int(b[i : i + 2], 16) for i in (0, 2, 4)]
    t = min(max(t, 0.0), 1.0)
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(ca, cb))


def window_color(w: Dict[str, Any]) -> str:
    """The genome overview's colour of one window.

    Trio: the share of origin-informative alleles from each parent, pink for
    the mother, blue for the father, grey when balanced. Duo: the missing
    parent's colour where the child lacks the present one's alleles. Solo,
    and a duo elsewhere: the het rate, dark where it falls to 0.
    """
    style = WINDOW_STYLE
    no_data = style.get("no_data", "#eeeeee")
    mat_side = [w.get("n_mat"), w.get("n_lacks_father")]
    pat_side = [w.get("n_pat"), w.get("n_lacks_mother")]
    has_mat = all(v is not None for v in mat_side)
    has_pat = all(v is not None for v in pat_side)
    minimum = int(style.get("balance_min_sites", 20))
    if has_mat and has_pat:
        mat, pat = sum(mat_side), sum(pat_side)
        if mat + pat < minimum:
            return no_data
        share = mat / (mat + pat)
        if share >= 0.5:
            return _mix(
                "#e5e7eb", style.get("maternal_error", "#ad1457"), (share - 0.5) * 2
            )
        return _mix(
            "#e5e7eb", style.get("paternal_error", "#0d47a1"), (0.5 - share) * 2
        )
    sites = w.get("n_sites") or 0
    if has_pat and sites and (w.get("n_lacks_mother") or 0) / sites >= 0.05:
        return style.get("paternal_error", "#0d47a1")
    if has_mat and sites and (w.get("n_lacks_father") or 0) / sites >= 0.05:
        return style.get("maternal_error", "#ad1457")
    rate = _het_rate(w)
    if rate is None:
        return no_data
    # 0.5 or more is ordinary; 0 is a run without any het call.
    return _mix(
        style.get("solo_het_dark", "#4a148c"),
        style.get("solo_het_light", "#f3e5f5"),
        rate / 0.5,
    )


def _runs(windows: Sequence[Dict[str, Any]]) -> List[Tuple[int, int, str]]:
    """Consecutive windows of one colour merged: ``(start, end, colour)``."""
    out: List[Tuple[int, int, str]] = []
    for w in windows:
        color = window_color(w)
        start, end = int(w["start"]), int(w["end"])
        if out and out[-1][2] == color and start <= out[-1][1] + 1:
            out[-1] = (out[-1][0], end, color)
        else:
            out.append((start, end, color))
    return out


def genome_svg(
    lanes: Sequence[str],
    windows: Dict[str, Dict[str, List[Dict[str, Any]]]],
    events: Dict[str, List[Dict[str, Any]]],
    modes: Dict[str, str],
) -> str:
    """Chromosomes 1-22 and X, one lane per member, its windows coloured.

    Args:
        lanes: Member IDs, in lane order.
        windows: Per member, per chromosome, its windows in position order.
        events: Per member, its events (outlined by confidence).
        modes: Per member, its analysis mode, shown in the lane's tooltip.
    """
    svg_w, lbl_w, right = 1200, 40, 25
    chrom_h, lane_h, lane_gap, row_gap = 8, 7, 1, 8
    max_mb = max(chrom_length_mb(c) for c in CHROMOSOMES)
    px_per_mb = (svg_w - lbl_w - right) / max_mb
    row_h = chrom_h + 2 + len(lanes) * (lane_h + lane_gap) + row_gap
    axis_y = len(CHROMOSOMES) * row_h
    legend_y = axis_y + 34
    svg_h = legend_y + 20
    no_data = WINDOW_STYLE.get("no_data", "#eeeeee")

    parts = [
        f'<svg viewBox="0 0 {svg_w} {svg_h}" xmlns="http://www.w3.org/2000/svg" '
        f'preserveAspectRatio="xMinYMin meet" '
        f'style="font-family: sans-serif; width: 100%; height: auto;">'
    ]
    parts += mb_axis_svg(lbl_w, px_per_mb, 0, axis_y, max_mb)
    for ci, chrom in enumerate(CHROMOSOMES):
        top = ci * row_h
        length_px = chrom_length_mb(chrom) * px_per_mb
        parts.append(
            f'<text x="{lbl_w - 6}" y="{top + chrom_h}" text-anchor="end" '
            f'font-size="11" fill="#374151">{chrom}</text>'
        )
        parts += chromosome_svg(chrom, lbl_w, top, px_per_mb, chrom_h)
        for li, iid in enumerate(lanes):
            y = top + chrom_h + 2 + li * (lane_h + lane_gap)
            mode = modes.get(iid, "")
            parts.append(
                f'<rect x="{lbl_w}" y="{y}" width="{length_px:.1f}" '
                f'height="{lane_h}" fill="{no_data}">'
                f"<title>{html.escape(iid)} ({html.escape(mode)})</title></rect>"
            )
            for start, end, color in _runs(windows.get(iid, {}).get(chrom, [])):
                x = lbl_w + start / 1e6 * px_per_mb
                w = max((end - start) / 1e6 * px_per_mb, 0.6)
                parts.append(
                    f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" '
                    f'height="{lane_h}" fill="{color}"/>'
                )
            for e in events.get(iid, []):
                if e.get("chrom") != chrom or e.get("start") is None:
                    continue
                x = lbl_w + int(e["start"]) / 1e6 * px_per_mb
                w = max((int(e["end"]) - int(e["start"])) / 1e6 * px_per_mb, 1.0)
                info = confidence_info(e.get("confidence"))
                tip = (
                    f"{iid} · chr{chrom} {span_text(e['start'], e['end'])} · "
                    f"{origin_info(e.get('origin'))['label']} "
                    f"{composition_label(e.get('composition'))} · {info['label']}"
                )
                parts.append(
                    f'<rect x="{x:.1f}" y="{y - 0.5}" width="{w:.1f}" '
                    f'height="{lane_h + 1}" fill="none" stroke="{info["color"]}" '
                    f'stroke-width="1.4"><title>{html.escape(tip)}</title></rect>'
                )
    parts.append(
        f'<text x="{svg_w / 2}" y="{axis_y + 28}" text-anchor="middle" '
        f'font-size="12" fill="#6b7280">Position (Mb)</text>'
    )
    legend = [
        ("Alleles from the mother only", WINDOW_STYLE.get("maternal_error", "#ad1457")),
        ("from the father only", WINDOW_STYLE.get("paternal_error", "#0d47a1")),
        ("balanced", "#e5e7eb"),
        ("no het call (solo)", WINDOW_STYLE.get("solo_het_dark", "#4a148c")),
        ("no data", no_data),
    ]
    lx = lbl_w
    for label, color in legend:
        parts.append(
            f'<rect x="{lx}" y="{legend_y}" width="14" height="10" rx="2" '
            f'fill="{color}" stroke="#d1d5db" stroke-width="0.5"/>'
        )
        parts.append(
            f'<text x="{lx + 18}" y="{legend_y + 9}" font-size="11" '
            f'fill="#374151">{html.escape(label)}</text>'
        )
        lx += 18 + len(label) * 6.2 + 18
    for key in ("high", "medium", "low"):
        info = confidence_info(key)
        parts.append(
            f'<rect x="{lx}" y="{legend_y}" width="14" height="10" fill="none" '
            f'stroke="{info["color"]}" stroke-width="1.4"/>'
        )
        text = f"{info['label']} event"
        parts.append(
            f'<text x="{lx + 18}" y="{legend_y + 9}" font-size="11" '
            f'fill="#374151">{html.escape(text)}</text>'
        )
        lx += 18 + len(text) * 6.2 + 18
    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Cohort
# ---------------------------------------------------------------------------


def cohort_chromosome_option(events: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """Events per chromosome, stacked by confidence."""
    counts: Dict[str, Dict[str, int]] = {}
    for e in events:
        counts.setdefault(e.get("confidence") or "", {}).setdefault(
            e.get("chrom") or "", 0
        )
        counts[e.get("confidence") or ""][e.get("chrom") or ""] += 1
    order = [c for c in CONFIDENCES if c in counts] + [
        c for c in counts if c not in CONFIDENCES
    ]
    series = [
        {
            "type": "bar",
            "name": confidence_info(conf)["label"],
            "stack": "events",
            "itemStyle": {"color": confidence_info(conf)["color"]},
            "data": [counts[conf].get(c, 0) for c in CHROMOSOMES],
        }
        for conf in order
    ]
    return {
        "animation": False,
        "legend": {"top": 0, "textStyle": {"fontSize": 11}},
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "shadow"}},
        "grid": {"left": 40, "right": 12, "top": 30, "bottom": 30},
        "xAxis": {
            "type": "category",
            "data": CHROMOSOMES,
            "axisLabel": {"fontSize": 10, "interval": 0},
        },
        "yAxis": {"type": "value", "minInterval": 1, "name": "events"},
        "series": series,
    }
