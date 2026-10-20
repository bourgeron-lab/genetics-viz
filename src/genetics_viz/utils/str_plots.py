"""ECharts builders for the STR views.

Pure functions returning option dicts, shared by the family tab and the cohort
panel so neither imports the other's privates.
"""

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .str_eh import status_info

#: (threshold column, label, class whose colour draws it)
_THRESHOLDS = (
    ("normal_max", "normal ≤", "normal"),
    ("intermediate_min", "intermediate ≥", "intermediate"),
    ("premutation_min", "premutation ≥", "premutation"),
    ("pathogenic_min", "pathogenic ≥", "pathogenic"),
)

#: Label position per threshold line. The intermediate threshold often sits
#: one unit past normal (HTT: 26 / 27); written along its line rather than
#: above it, its label does not overprint the normal one.
_LABEL_POSITION = {"intermediate_min": "insideEndTop"}


def _num(value: Any) -> Optional[float]:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def allele_class(count: float, thresholds: Dict[str, Any]) -> str:
    """Class of a repeat count, in the pipeline's order (eh_family_table)."""
    normal = _num(thresholds.get("normal_max"))
    inter = _num(thresholds.get("intermediate_min"))
    premut = _num(thresholds.get("premutation_min"))
    patho = _num(thresholds.get("pathogenic_min"))
    if normal is None and inter is None and premut is None and patho is None:
        return "unclassified"
    if patho is not None and count >= patho:
        return "pathogenic"
    if premut is not None and count >= premut:
        return "premutation"
    if inter is not None and count >= inter:
        return "intermediate"
    if normal is None or count <= normal:
        return "normal"
    return "above_normal"


def distribution_option(
    counts: Dict[int, int],
    thresholds: Dict[str, Any],
    members: Sequence[Tuple[str, Sequence[int], str]],
    title: str = "",
) -> Dict[str, Any]:
    """Histogram of the founders' longest allele, members marked above it.

    Args:
        counts: repeat count -> number of founders (``locus_distribution``).
        thresholds: a row holding ``normal_max``, ``intermediate_min``,
            ``premutation_min``, ``pathogenic_min``.
        members: ``(label, alleles, status)`` for each member to mark.
        title: chart title.

    The x-axis stops a little past the longest allele seen, so a threshold
    far out (YEATS2's pathogenic 1000) does not crush the histogram: it is
    named in the subtitle instead. The count axis is logarithmic, since one
    common allele holds most founders and the tail is what matters.
    """
    seen = list(counts) + [a for _, alleles, _ in members for a in alleles]
    hi = max(seen) if seen else 10
    lo = min(seen) if seen else 0
    limits = {
        col: _num(thresholds.get(col))
        for col, _, _ in _THRESHOLDS
        if _num(thresholds.get(col)) is not None
    }
    x_max = math.ceil(max(hi * 1.15, hi + 3))
    near = [v for v in limits.values() if v <= x_max * 1.5]
    if near:
        x_max = math.ceil(max(x_max, max(near) * 1.1))
    x_min = max(0, math.floor(min(lo, *near) if near else lo) - 2)

    off_axis = [
        f"{label} {v:g}"
        for col, label, _ in _THRESHOLDS
        if (v := limits.get(col)) is not None and v > x_max
    ]
    subtitle = "Bars: the cohort founders' longest allele (log scale)."
    if off_axis:
        subtitle += " Off the axis: " + ", ".join(off_axis) + "."

    bars = [
        {
            "value": [x, n],
            "itemStyle": {"color": status_info(allele_class(x, thresholds))["color"]},
        }
        for x, n in sorted(counts.items())
    ]
    mark_lines = [
        {
            "xAxis": v,
            "name": f"{label} {v:g}",
            "lineStyle": {"color": status_info(cls)["color"], "type": "dashed"},
            "label": {
                "formatter": f"{label} {v:g}",
                "position": _LABEL_POSITION.get(col, "end"),
            },
        }
        for col, label, cls in _THRESHOLDS
        if (v := limits.get(col)) is not None and v <= x_max
    ]

    labels = [m[0] for m in members]
    points: List[Dict[str, Any]] = []
    for label, alleles, _status in members:
        for a in alleles:
            points.append(
                {
                    "value": [a, label],
                    "itemStyle": {
                        "color": status_info(allele_class(a, thresholds))["color"]
                    },
                }
            )
    member_height = max(60, 22 * len(labels))

    return {
        "title": {
            "text": title,
            "subtext": subtitle,
            "left": "center",
            "textStyle": {"fontSize": 14},
            "subtextStyle": {"fontSize": 11},
        },
        "tooltip": {"trigger": "item"},
        "grid": [
            {"left": 140, "right": 40, "top": 60, "height": member_height},
            {"left": 140, "right": 40, "top": 80 + member_height, "bottom": 50},
        ],
        "xAxis": [
            {
                "type": "value",
                "gridIndex": 0,
                "min": x_min,
                "max": x_max,
                "axisLabel": {"show": False},
                "splitLine": {"show": True, "lineStyle": {"color": "#eee"}},
            },
            {
                "type": "value",
                "gridIndex": 1,
                "min": x_min,
                "max": x_max,
                "name": "Repeat units",
                "nameLocation": "middle",
                "nameGap": 28,
            },
        ],
        "yAxis": [
            {
                "type": "category",
                "gridIndex": 0,
                "data": labels,
                "inverse": True,
                "axisTick": {"show": False},
            },
            {
                "type": "log",
                "gridIndex": 1,
                "min": 1,
                "name": "Founders",
                "nameLocation": "middle",
                "nameGap": 40,
            },
        ],
        "series": [
            {
                "type": "scatter",
                "name": "Members",
                "xAxisIndex": 0,
                "yAxisIndex": 0,
                "symbol": "diamond",
                "symbolSize": 14,
                "data": points,
                "tooltip": {"formatter": "{c}"},
            },
            {
                "type": "bar",
                "name": "Founders",
                "xAxisIndex": 1,
                "yAxisIndex": 1,
                "barMaxWidth": 12,
                "data": bars,
                "markLine": {"symbol": "none", "data": mark_lines},
                "tooltip": {"formatter": "{c}"},
            },
        ],
    }
