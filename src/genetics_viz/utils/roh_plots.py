"""Shared plot builders for the ROH tabs.

Pure functions returning ECharts option dicts or SVG strings - nothing here
touches ``ui`` - used by both the per-family ROH tab and the cohort one. Every
plot puts a person next to the unrelated reference genomes of their own
ancestry region; the fixed thresholds come second, as lines.
"""

import html
import math
import random
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .cytobands import CHROM_ORDER
from .ideogram import chrom_length_mb, chromosome_svg, mb_axis_svg
from .roh import (
    CLASSES,
    EXCESS_PERCENTILE,
    class_color,
    pattern_info,
    threshold_grades,
)

#: Values are rounded before serialisation; the reference adds ~750 points
#: per region to each plot.
_PRECISION = 5

_REFERENCE_COLOR = "rgba(120, 120, 120, 0.35)"

#: The ROH panel covers the autosomes only.
AUTOSOMES = [c for c in CHROM_ORDER if c.isdigit()]


def _round(value: Optional[float]) -> Optional[float]:
    return None if value is None else round(value, _PRECISION)


def _axis_max(values: Sequence[float], lines: Sequence[float]) -> float:
    """Upper bound of a value axis: the data, then the next threshold above it.

    Showing the next fixed threshold gives scale to the data without squeezing
    a European distribution (max ~0.09) into the first tenth of a 0-0.17 axis.
    """
    top = max(values, default=0.0)
    above = [t for t in lines if t > top]
    bound = min(above) if above else top
    # Round up to a whole hundredth so the end tick reads 0.19, not 0.1985.
    return max(math.ceil(bound * 1.08 * 100) / 100, 0.01)


def threshold_marklines(show_relationship: bool, axis: str = "x") -> List[Dict]:
    """Dashed lines at the fixed F_ROH thresholds, labelled with the grade."""
    lines = []
    for grade in threshold_grades():
        name = grade.grade
        if show_relationship:
            name = f"{grade.grade} ({grade.relationship.split(':')[0]})"
        lines.append(
            {
                axis + "Axis": grade.threshold,
                "name": name,
                "lineStyle": {"color": grade.color, "type": "dashed", "width": 1.5},
                "label": {
                    "formatter": f"{name} > {grade.threshold:g}",
                    "color": grade.color,
                    "fontSize": 10,
                    # Written along the line: three labels side by side above
                    # the plot collide, and the curator wording is long.
                    "position": "insideEndTop",
                },
            }
        )
    return lines


def _percentile_markline(value: float, axis: str, label: str) -> Dict:
    return {
        axis + "Axis": _round(value),
        "name": label,
        "lineStyle": {"color": "#455a64", "type": "dotted", "width": 1.5},
        "label": {
            "formatter": label,
            "color": "#455a64",
            "fontSize": 10,
            # A vertical line's label goes at the bottom: the threshold labels
            # sit at the top, and the 95th percentile is often close to one.
            "position": "insideEndTop" if axis == "y" else "insideStartTop",
        },
    }


def reference_strip_option(
    reference: List[float],
    members: List[Dict[str, Any]],
    *,
    title: str,
    metric_label: str,
    p95: Optional[float],
    show_thresholds: bool,
    show_relationship: bool,
) -> Dict:
    """Horizontal strip: the reference region jittered, the members as diamonds.

    Args:
        reference: Reference values of the comparison group.
        members: ``{"name", "value", "color"}`` per member to highlight.
        title: Plot title, naming the comparison group and its size.
        metric_label: Axis name.
        p95: The group's 95th percentile, drawn dotted; ``None`` to omit.
        show_thresholds: Draw the fixed F_ROH thresholds (only meaningful for
            ``froh_1.5mb``).
        show_relationship: Curator wording on the threshold labels.
    """
    rng = random.Random(0)
    ref_data = [[_round(v), round(rng.uniform(-1, 1), 3)] for v in reference]

    n = len(members)
    member_data = []
    for i, m in enumerate(members):
        # Spread members vertically so several at similar values stay apart.
        y = 0.0 if n == 1 else -0.6 + 1.2 * i / (n - 1)
        member_data.append(
            {
                "name": m["name"],
                "value": [_round(m["value"]), round(y, 3)],
                "itemStyle": {
                    "color": m.get("color", "#1f2937"),
                    "borderColor": "#1f2937",
                    "borderWidth": 1.5,
                },
            }
        )

    mark_lines: List[Dict] = []
    thresholds: List[float] = []
    if show_thresholds:
        mark_lines += threshold_marklines(show_relationship, axis="x")
        thresholds = [g.threshold for g in threshold_grades()]
    if p95 is not None:
        mark_lines.append(
            _percentile_markline(p95, "x", f"{EXCESS_PERCENTILE:g}th pct")
        )

    values = list(reference) + [m["value"] for m in members]
    return {
        "title": {
            "text": title,
            "left": "center",
            "textStyle": {"fontSize": 12, "fontWeight": "normal"},
        },
        "tooltip": {
            "trigger": "item",
            ":formatter": (
                "p => p.seriesName === 'Reference' ? "
                "'Reference: ' + p.value[0].toFixed(4) : "
                "p.name + ': ' + p.value[0].toFixed(4)"
            ),
        },
        # Room above the plot for the title and, under it, the threshold labels.
        "grid": {"top": 50, "bottom": 40, "left": 20, "right": 30},
        "xAxis": {
            "type": "value",
            "name": metric_label,
            "nameLocation": "middle",
            "nameGap": 25,
            "min": 0,
            # The y range spans zero, so the axis line would otherwise be drawn
            # through the middle of the strip.
            "axisLine": {"onZero": False},
            "max": _axis_max(values, thresholds),
        },
        "yAxis": {"type": "value", "show": False, "min": -1.3, "max": 1.3},
        "series": [
            {
                "name": "Reference",
                "type": "scatter",
                "data": ref_data,
                "symbolSize": 4,
                "itemStyle": {"color": _REFERENCE_COLOR},
                "large": True,
                "markLine": {"symbol": "none", "silent": True, "data": mark_lines},
            },
            {
                "name": "Members",
                "type": "scatter",
                "data": member_data,
                "symbol": "diamond",
                "symbolSize": 14,
                "z": 10,
                "label": {
                    "show": True,
                    "position": "top",
                    "fontSize": 10,
                    "formatter": "{b}",
                },
                "labelLayout": {"hideOverlap": True},
            },
        ],
    }


def pattern_scatter_option(
    reference: List[Tuple[float, float]],
    points: List[Dict[str, Any]],
    *,
    title: str,
    x_line: Optional[float],
    y_line: Optional[float],
    show_labels: bool = True,
    symbol_size: int = 14,
) -> Dict:
    """Population background against recent: the ``homozygosity_pattern`` plane.

    x is ``f_short + f_intermediate`` (old shared ancestry), y is ``f_long``
    (a recent common ancestor of the parents). The dotted lines are the
    comparison group's 95th percentiles; the four quadrants they cut are the
    four patterns the pipeline reports.

    Args:
        reference: ``(x, y)`` of every reference individual of the group.
        points: ``{"name", "x", "y", "color"}`` per highlighted individual.
        x_line: The group's 95th percentile of the background.
        y_line: The threshold for "recent": the 95th percentile of ``f_long``,
            floored at the tool's minimum.
    """
    mark_lines = []
    if x_line is not None:
        mark_lines.append(
            _percentile_markline(x_line, "x", f"{EXCESS_PERCENTILE:g}th pct")
        )
    if y_line is not None:
        mark_lines.append(
            _percentile_markline(y_line, "y", f"{EXCESS_PERCENTILE:g}th pct")
        )

    point_data = [
        {
            "name": p["name"],
            "value": [_round(p["x"]), _round(p["y"])],
            "itemStyle": {
                "color": p.get("color", "#1f2937"),
                "borderColor": "#1f2937",
                "borderWidth": 1 if symbol_size >= 10 else 0.5,
            },
        }
        for p in points
    ]

    corner = {"fontSize": 11, "fontWeight": "bold"}

    def quadrant(text: str, color: str, **pos: Any) -> Dict:
        return {
            "type": "text",
            "silent": True,
            "style": {"text": text, "fill": color, **corner},
            **pos,
        }

    return {
        "title": {
            "text": title,
            "left": "center",
            "textStyle": {"fontSize": 12, "fontWeight": "normal"},
        },
        "tooltip": {
            "trigger": "item",
            ":formatter": (
                "p => (p.seriesName === 'Reference' ? 'Reference' : p.name)"
                " + '<br>short+intermediate: ' + p.value[0].toFixed(4)"
                " + '<br>long: ' + p.value[1].toFixed(4)"
            ),
        },
        "grid": {"top": 50, "bottom": 50, "left": 60, "right": 30},
        "graphic": [
            quadrant(
                pattern_info("recent")["label"],
                pattern_info("recent")["color"],
                left=70,
                top=45,
            ),
            quadrant(
                pattern_info("recent+population")["label"],
                pattern_info("recent+population")["color"],
                right=40,
                top=45,
            ),
            quadrant(
                pattern_info("population")["label"],
                pattern_info("population")["color"],
                right=40,
                bottom=60,
            ),
        ],
        "xAxis": {
            "type": "value",
            "name": "F short + intermediate (population background)",
            "nameLocation": "middle",
            "nameGap": 28,
            "min": 0,
        },
        "yAxis": {
            "type": "value",
            "name": "F long (recent)",
            "nameLocation": "middle",
            "nameGap": 45,
            "min": 0,
        },
        "series": [
            {
                "name": "Reference",
                "type": "scatter",
                "data": [[_round(x), _round(y)] for x, y in reference],
                "symbolSize": 4,
                "itemStyle": {"color": _REFERENCE_COLOR},
                "large": True,
                "markLine": {"symbol": "none", "silent": True, "data": mark_lines},
            },
            {
                "name": "Individuals",
                "type": "scatter",
                "data": point_data,
                "symbol": "diamond",
                "symbolSize": symbol_size,
                "z": 10,
                "label": {
                    "show": show_labels,
                    "position": "top",
                    "fontSize": 10,
                    "formatter": "{b}",
                },
                "labelLayout": {"hideOverlap": True},
            },
        ],
    }


# ---------------------------------------------------------------------------
# Genome view
# ---------------------------------------------------------------------------

#: One segment: ``(chrom, start_bp, end_bp, roh_class, tooltip)``.
Segment = Tuple[str, int, int, str, str]
#: One hotspot interval: ``(chrom, start_bp, end_bp)``.
Interval = Tuple[str, int, int]


def genome_svg(
    lanes: List[str],
    segments: Dict[str, List[Segment]],
    hotspots: Dict[str, List[Interval]],
    upd: Dict[str, str],
) -> str:
    """Chromosomes 1-22, each with one thin lane per individual underneath.

    Args:
        lanes: Individual IDs, in lane order.
        segments: Per individual, their ROH segments.
        hotspots: Per individual, the reference hotspots of their region,
            shaded in their lane so a segment sitting in one stands out as
            less informative.
        upd: Per individual, a UPD candidate chromosome, outlined in red.
    """
    svg_w, lbl_w, right = 1200, 40, 25
    chrom_h, lane_h, lane_gap, row_gap = 8, 6, 1, 8
    max_mb = max(chrom_length_mb(c) for c in AUTOSOMES)
    px_per_mb = (svg_w - lbl_w - right) / max_mb
    row_h = chrom_h + 2 + len(lanes) * (lane_h + lane_gap) + row_gap
    axis_y = len(AUTOSOMES) * row_h
    legend_y = axis_y + 34
    svg_h = legend_y + 20

    by_chrom: Dict[Tuple[str, str], List[Segment]] = {}
    for iid, segs in segments.items():
        for seg in segs:
            by_chrom.setdefault((iid, seg[0]), []).append(seg)
    hot_by_chrom: Dict[Tuple[str, str], List[Interval]] = {}
    for iid, intervals in hotspots.items():
        for iv in intervals:
            hot_by_chrom.setdefault((iid, iv[0]), []).append(iv)

    parts = [
        f'<svg viewBox="0 0 {svg_w} {svg_h}" xmlns="http://www.w3.org/2000/svg" '
        f'preserveAspectRatio="xMinYMin meet" '
        f'style="font-family: sans-serif; width: 100%; height: auto;">'
    ]
    parts += mb_axis_svg(lbl_w, px_per_mb, 0, axis_y, max_mb)

    for ci, chrom in enumerate(AUTOSOMES):
        top = ci * row_h
        length_px = chrom_length_mb(chrom) * px_per_mb
        parts.append(
            f'<text x="{lbl_w - 6}" y="{top + chrom_h}" text-anchor="end" '
            f'font-size="11" fill="#374151">{chrom}</text>'
        )
        parts += chromosome_svg(chrom, lbl_w, top, px_per_mb, chrom_h)
        for li, iid in enumerate(lanes):
            y = top + chrom_h + 2 + li * (lane_h + lane_gap)
            outline = (
                ' stroke="#d32f2f" stroke-width="1.2"' if upd.get(iid) == chrom else ""
            )
            parts.append(
                f'<rect x="{lbl_w}" y="{y}" width="{length_px:.1f}" '
                f'height="{lane_h}" fill="#f8fafc"{outline}>'
                f"<title>{html.escape(iid)}</title></rect>"
            )
            for _, start, end in hot_by_chrom.get((iid, chrom), []):
                x = lbl_w + start / 1e6 * px_per_mb
                w = max((end - start) / 1e6 * px_per_mb, 0.8)
                parts.append(
                    f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{lane_h}" '
                    f'fill="#cbd5e1"><title>Reference ROH hotspot</title></rect>'
                )
            for _, start, end, cls, tip in by_chrom.get((iid, chrom), []):
                x = lbl_w + start / 1e6 * px_per_mb
                w = max((end - start) / 1e6 * px_per_mb, 0.8)
                parts.append(
                    f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="{lane_h}" '
                    f'fill="{class_color(cls)}"><title>{html.escape(tip)}</title></rect>'
                )

    parts.append(
        f'<text x="{svg_w / 2}" y="{axis_y + 28}" text-anchor="middle" '
        f'font-size="12" fill="#6b7280">Position (Mb)</text>'
    )
    lx = lbl_w
    legend = [(info["label"], info["color"]) for info in CLASSES.values()]
    legend.append(("Reference hotspot", "#cbd5e1"))
    for label, color in legend:
        parts.append(
            f'<rect x="{lx}" y="{legend_y}" width="14" height="10" rx="2" fill="{color}"/>'
        )
        parts.append(
            f'<text x="{lx + 18}" y="{legend_y + 9}" font-size="11" '
            f'fill="#374151">{html.escape(label)}</text>'
        )
        lx += 18 + len(label) * 6.5 + 18
    parts.append(
        f'<rect x="{lx}" y="{legend_y}" width="14" height="10" fill="#f8fafc" '
        f'stroke="#d32f2f" stroke-width="1.2"/>'
    )
    parts.append(
        f'<text x="{lx + 18}" y="{legend_y + 9}" font-size="11" '
        f'fill="#374151">UPD candidate chromosome</text>'
    )
    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Cohort ROH islands
# ---------------------------------------------------------------------------


def _offsets() -> Dict[str, float]:
    """Cumulative start of each autosome on a genome-wide Mb axis."""
    out, total = {}, 0.0
    for chrom in AUTOSOMES:
        out[chrom] = total
        total += chrom_length_mb(chrom)
    return out


def islands_option(
    windows: List[Tuple[str, int, int, float]],
    hotspots: List[Interval],
    *,
    title: str,
) -> Dict:
    """Genome-wide frequency of founders in ROH, with the reference hotspots.

    Args:
        windows: ``(chrom, start, end, freq)`` per 100 kb window. Only windows
            with at least one founder in ROH are written, so gaps are zeros.
        hotspots: The reference ROH islands of the same group, shaded.
    """
    offsets = _offsets()
    data: List[List[Any]] = []
    prev: Optional[Tuple[str, int]] = None
    for chrom, start, end, freq in sorted(
        (w for w in windows if w[0] in offsets),
        key=lambda w: (AUTOSOMES.index(w[0]), w[1]),
    ):
        x0 = offsets[chrom] + start / 1e6
        # Drop to zero across a gap, rather than drawing a slope over it.
        if prev is None or prev != (chrom, start):
            if prev is not None:
                pchrom, pend = prev
                data.append(
                    [
                        round(offsets[pchrom] + pend / 1e6, 2),
                        0,
                        pchrom,
                        round(pend / 1e6, 1),
                    ]
                )
            data.append([round(x0, 2), 0, chrom, round(start / 1e6, 1)])
        data.append([round(x0 + 0.05, 2), freq, chrom, round(start / 1e6, 1)])
        prev = (chrom, end)

    top = max((row[1] for row in data), default=0.0)
    marker_y = round(top * 1.06 + 0.005, 3)
    hot_areas = [
        [
            {"xAxis": round(offsets[c] + s / 1e6, 2)},
            {"xAxis": round(offsets[c] + e / 1e6, 2)},
        ]
        for c, s, e in hotspots
        if c in offsets
    ]
    # Alternate chromosome shading, labelled, on an empty helper series.
    chrom_areas = [
        [
            {
                "xAxis": round(offsets[c], 2),
                "itemStyle": {
                    "color": "rgba(0,0,0,0.035)" if i % 2 else "rgba(0,0,0,0)"
                },
                "label": {
                    "show": True,
                    "formatter": c,
                    "position": "insideTop",
                    "fontSize": 10,
                    "color": "#6b7280",
                },
            },
            {"xAxis": round(offsets[c] + chrom_length_mb(c), 2)},
        ]
        for i, c in enumerate(AUTOSOMES)
    ]
    genome_mb = round(sum(chrom_length_mb(c) for c in AUTOSOMES), 1)
    return {
        "title": {
            "text": title,
            "left": "center",
            "textStyle": {"fontSize": 13, "fontWeight": "normal"},
        },
        "tooltip": {
            "trigger": "axis",
            ":formatter": (
                "ps => { const p = ps.find(q => q.seriesName === 'Founders in ROH');"
                " const h = ps.find(q => q.seriesName === 'Reference hotspots');"
                " if (!p) return h ? 'Reference hotspot, chr' + h.value[2] + ':'"
                " + h.value[3] + ' Mb' : '';"
                " return 'chr' + p.value[2] + ':' + p.value[3] + ' Mb<br>'"
                " + (100 * p.value[1]).toFixed(1) + '% of founders'"
                " + (h ? '<br>near a reference hotspot' : ''); }"
            ),
        },
        "legend": {"bottom": 0, "data": ["Founders in ROH", "Reference hotspots"]},
        "grid": {"top": 50, "bottom": 70, "left": 60, "right": 20},
        "dataZoom": [
            {"type": "inside"},
            {"type": "slider", "bottom": 28, "height": 16},
        ],
        "xAxis": {
            "type": "value",
            "min": 0,
            "max": genome_mb,
            "axisLabel": {"show": False},
            "splitLine": {"show": False},
        },
        "yAxis": {
            "type": "value",
            "name": "Fraction of founders in ROH",
            "nameLocation": "middle",
            "nameGap": 45,
            "min": 0,
        },
        "series": [
            {
                "name": "Chromosomes",
                "type": "line",
                "data": [],
                "silent": True,
                "markArea": {"silent": True, "data": chrom_areas},
            },
            {
                "name": "Founders in ROH",
                "type": "line",
                "data": data,
                "showSymbol": False,
                "sampling": "lttb",
                "lineStyle": {"width": 1, "color": "#1e88e5"},
                "areaStyle": {"color": "rgba(30,136,229,0.2)"},
                "itemStyle": {"color": "#1e88e5"},
            },
            {
                # A hotspot is 0.1-1 Mb of a 2,875 Mb axis: under a pixel at
                # full zoom. A marker row above the curve shows where they are;
                # the shaded areas take over once zoomed in.
                "name": "Reference hotspots",
                "type": "scatter",
                "data": [
                    [
                        round(offsets[c] + (s + e) / 2e6, 2),
                        marker_y,
                        c,
                        round(s / 1e6, 1),
                    ]
                    for c, s, e in hotspots
                    if c in offsets
                ],
                "symbol": "triangle",
                "symbolRotate": 180,
                "symbolSize": 7,
                "itemStyle": {"color": "rgba(216,27,96,0.8)"},
                "markArea": {
                    "silent": True,
                    "itemStyle": {"color": "rgba(216,27,96,0.18)"},
                    "data": hot_areas,
                },
            },
        ],
    }
