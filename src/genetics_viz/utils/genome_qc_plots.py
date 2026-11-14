"""ECharts builder for the genome QC views (family and cohort Sample QC tabs).

One metric per chart, one row per ``seq_class``: the pipeline judges a sample
against the samples of its own class, so the chart shows it that way. Pure
functions returning option dicts.
"""

import html
import json
import math
import random
import statistics
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .genome_qc import SEQ_CLASSES, THRESHOLDS, MetricInfo, seq_class_info
from .roh import to_float

_BACKGROUND = "#b0bec5"
_FLAGGED = "#c62828"
_MEMBER = "#1565c0"
_LINE = "#455a64"
_TITLE_STYLE = {"fontSize": 12, "fontWeight": "normal"}

#: A log axis cannot draw 0: an EBV count of zero is drawn here.
LOG_FLOOR = 1e-4

#: Scaled MAD, as genomeqc_tables computes it.
MAD_SCALE = 1.4826


def _esc(text: Any) -> str:
    return html.escape(str(text), quote=True)


def class_order(classes: Iterable[Optional[str]]) -> List[str]:
    """The classes present, in ``genome_qc.yaml`` order, then any other."""
    present = {c or "unknown" for c in classes}
    known = [c for c in SEQ_CLASSES if c in present]
    return known + sorted(present - set(known))


def class_fences(
    rows: Sequence[Dict[str, Any]], metric: MetricInfo
) -> Dict[str, Tuple[float, List[float]]]:
    """Per class: the median and the outlier fence(s) on the metric's side.

    The same rule as ``genomeqc_tables`` (median ± 5 scaled MADs, 20 samples
    at least), drawn for orientation; the flags are the pipeline's.
    """
    if not metric.side:
        return {}
    k = THRESHOLDS.get("mad_k", 5.0)
    minimum = int(THRESHOLDS.get("min_class", 20))
    by_class: Dict[str, List[float]] = {}
    for r in rows:
        value = to_float(r.get(metric.key))
        if value is not None:
            by_class.setdefault(r.get("seq_class") or "unknown", []).append(value)
    out = {}
    for cls, values in by_class.items():
        if len(values) < minimum:
            continue
        med = statistics.median(values)
        mad = MAD_SCALE * statistics.median(abs(v - med) for v in values)
        if mad <= 0:
            continue
        fences = []
        if metric.side in ("low", "both"):
            fences.append(med - k * mad)
        if metric.side in ("high", "both"):
            fences.append(med + k * mad)
        out[cls] = (med, fences)
    return out


def class_strip_option(
    metric: MetricInfo,
    rows: Sequence[Dict[str, Any]],
    *,
    highlighted: Set[str] = frozenset(),
    flagged: Set[str] = frozenset(),
    reference: Optional[Sequence[Dict[str, Any]]] = None,
) -> Dict:
    """One metric, one jittered strip per ``seq_class``.

    Args:
        rows: Samples-table rows to draw (``sample``, ``family_id``,
            ``seq_class`` and the metric's column).
        highlighted: Samples drawn as labelled diamonds (a family's members);
            the rest are small points.
        flagged: Samples the pipeline flagged on this metric, drawn red.
        reference: Rows the class fences are computed from (the whole
            cohort); default ``rows``.
    """
    classes = class_order(r.get("seq_class") for r in rows)
    index = {c: i for i, c in enumerate(classes)}
    rng = random.Random(0)
    background, red, members = [], [], []
    for r in rows:
        value = to_float(r.get(metric.key))
        if value is None:
            continue
        sample = r.get("sample") or ""
        cls = r.get("seq_class") or "unknown"
        shown = max(value, LOG_FLOOR) if metric.log else value
        is_member = sample in highlighted
        y = index[cls] + rng.uniform(-0.3, 0.3) * (0.5 if is_member else 1.0)
        tip = (
            f"<b>{_esc(sample)}</b> ({_esc(r.get('family_id') or '')})<br>"
            f"{_esc(seq_class_info(cls)['label'])}<br>"
            f"{_esc(metric.label)}: {_esc(metric.fmt(value))}"
        )
        point = {"name": sample, "value": [shown, round(y, 3)], "tip": tip}
        if is_member:
            point["itemStyle"] = {
                "color": _FLAGGED if sample in flagged else _MEMBER,
                "borderColor": "#1f2937",
                "borderWidth": 1,
            }
            members.append(point)
        elif sample in flagged:
            red.append(point)
        else:
            background.append(point)

    lines: List[Dict[str, Any]] = []
    for value in metric.lines:
        lines.append(
            {
                "xAxis": value,
                "name": f"{value:g}",
                "lineStyle": {"color": _LINE, "type": "dashed", "width": 1},
                "label": {"formatter": f"{value:g}", "color": _LINE, "fontSize": 10},
            }
        )
    for cls, (med, fences) in class_fences(reference or rows, metric).items():
        if cls not in index:
            continue
        i = index[cls]
        for value, style in [(med, "solid")] + [(f, "dotted") for f in fences]:
            if metric.log and value <= 0:
                continue
            lines.append(
                [
                    {
                        "coord": [value, i - 0.42],
                        "lineStyle": {"color": _LINE, "type": style, "width": 1},
                        "label": {"show": False},
                    },
                    {"coord": [value, i + 0.42]},
                ]
            )

    labels = json.dumps([seq_class_info(c)["label"] for c in classes])
    x_axis: Dict[str, Any] = {
        "type": "log" if metric.log else "value",
        "name": metric.label,
        "nameLocation": "middle",
        "nameGap": 25,
        "scale": not metric.log,
        "axisLine": {"onZero": False},
    }
    if metric.log:
        # The axis reaches past the highest threshold (CELL_LINE's 1), so a
        # cohort of blood samples still shows where a cell line would sit.
        top = max(
            [p["value"][0] for p in background + red + members] + metric.lines,
            default=1.0,
        )
        x_axis["min"] = LOG_FLOOR
        x_axis["max"] = 10 ** math.ceil(math.log10(max(top, LOG_FLOOR) * 1.5))
    return {
        "title": {"text": metric.label, "left": "center", "textStyle": _TITLE_STYLE},
        "tooltip": {
            "trigger": "item",
            "confine": True,
            ":formatter": "p => (p.data && p.data.tip) || ''",
        },
        "grid": {"top": 48, "bottom": 40, "left": 90, "right": 30},
        "xAxis": x_axis,
        "yAxis": {
            "type": "value",
            "min": -0.5,
            "max": len(classes) - 0.5,
            # Ticks every half step, labelled on the whole ones: each class
            # label sits in the middle of its band.
            "interval": 0.5,
            "splitLine": {"show": False},
            "axisTick": {"show": False},
            "axisLabel": {
                ":formatter": (
                    f"v => Number.isInteger(v) ? (({labels})[v] ?? '') : ''"
                ),
            },
        },
        "series": [
            {
                "name": "Cohort",
                "type": "scatter",
                "data": background,
                "symbolSize": 5,
                "itemStyle": {"color": _BACKGROUND},
                "markLine": {"symbol": "none", "silent": True, "data": lines},
            },
            {
                "name": "Flagged",
                "type": "scatter",
                "data": red,
                "symbolSize": 7,
                "itemStyle": {"color": _FLAGGED},
                "z": 5,
            },
            {
                "name": "Members",
                "type": "scatter",
                "data": members,
                "symbol": "diamond",
                "symbolSize": 13,
                "z": 10,
                "label": {
                    "show": len(members) <= 12,
                    "position": "top",
                    "fontSize": 9,
                    "formatter": "{b}",
                },
                "labelLayout": {"hideOverlap": True},
            },
        ],
    }


def chart_height(rows: Sequence[Dict[str, Any]]) -> int:
    """Pixels for a chart of these rows' classes."""
    return 110 + 70 * len(class_order(r.get("seq_class") for r in rows))
