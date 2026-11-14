"""ECharts options of the cohort De novo tab: pure functions over the rows of
``denovo.child_rows`` and ``denovo.flag_summary``.

The literature lines come from ``expected`` in ``config/denovo.yaml``.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from . import denovo

_SNV = "#1976d2"
_INDEL = "#ef6c00"
_LIT = "rgba(46, 125, 50, 0.12)"
_LIT_LINE = "#2e7d32"
_QC_LINE = "#c62828"
_CELL_LINE = "#880e4f"
_CELL_LINE_INDEL = "#d81b60"


def _label(c: Dict[str, Any]) -> str:
    return f"{c['child']} ({c['FID']})"


def rate_option(children: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Strict SNVs and indels per callable Gb per child, the literature band.

    Bars stacked, sorted by the total; the green band is the published SNV
    rate (2 × μ), the red dashed lines the cohort's QC range.
    """
    rows = sorted(
        (c for c in children if c.get("rate_per_gb") is not None),
        key=lambda c: c["rate_per_gb"],
    )
    band = denovo.EXPECTED.get("snv_rate_per_gb") or []
    # One QC range per class (genome QC): drawn only when the children shown
    # share one, as a line across children of other classes would mislead.
    ranges = {
        (c.get("qc_low"), c.get("qc_high"))
        for c in rows
        if c.get("qc_low") is not None and c.get("qc_high") is not None
    }
    low, high = next(iter(ranges)) if len(ranges) == 1 else (None, None)
    scale_rows = [c for c in rows if not c.get("cell_line")] or rows

    def bar(c: Dict[str, Any], value: Any, color: str) -> Any:
        # A cell line's bars in its own colour: its count is mostly somatic.
        if not c.get("cell_line"):
            return value
        return {"value": value, "itemStyle": {"color": color}}

    snv_series: Dict[str, Any] = {
        "name": "SNVs",
        "type": "bar",
        "stack": "rate",
        "color": _SNV,
        "data": [bar(c, c["snv_per_gb"], _CELL_LINE) for c in rows],
    }
    if len(band) == 2:
        snv_series["markArea"] = {
            "silent": True,
            "itemStyle": {"color": _LIT},
            "label": {"show": True, "position": "insideTopLeft", "color": _LIT_LINE},
            "data": [
                [
                    {"name": "Published SNV rate", "yAxis": band[0]},
                    {"yAxis": band[1]},
                ]
            ],
        }
    lines = []
    for name, value in (("QC low", low), ("QC high", high)):
        if value:
            lines.append({"name": name, "yAxis": round(value, 1)})
    indel_series: Dict[str, Any] = {
        "name": "Indels",
        "type": "bar",
        "stack": "rate",
        "color": _INDEL,
        "data": [bar(c, c["indel_per_gb"], _CELL_LINE_INDEL) for c in rows],
    }
    if lines:
        indel_series["markLine"] = {
            "silent": True,
            "symbol": "none",
            "lineStyle": {"color": _QC_LINE, "type": "dashed"},
            "label": {"formatter": "{b}: {c}", "position": "insideEndTop"},
            "data": lines,
        }
    return {
        "title": {
            "text": "Strict candidates per callable Gb",
            "subtext": "Autosomes, unmasked callable genome; green: published SNV"
            " rate (2 × μ); red: the QC range"
            + ("" if len(ranges) <= 1 else " (one per class: see the table)")
            + (
                "; purple: cell lines (their bars may run past the axis)"
                if any(c.get("cell_line") for c in rows)
                else ""
            ),
            "left": "center",
            "textStyle": {"fontSize": 14},
        },
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "shadow"}},
        "legend": {"top": 62},
        "grid": {"top": 95, "left": 50, "right": 20, "bottom": 90},
        "xAxis": {
            "type": "category",
            "data": [_label(c) for c in rows],
            "axisLabel": {"rotate": 60, "fontSize": 10},
        },
        # The QC lines may lie above every bar: the axis reaches them. Cell
        # lines (up to 650 per Gb against a median of 28) would flatten the
        # others: the axis is set on the rest, and their bars run past it.
        "yAxis": {
            "type": "value",
            "name": "per Gb",
            "max": round(
                max(
                    [
                        high or 0,
                        *(c["rate_per_gb"] for c in scale_rows),
                        *band,
                    ]
                )
                * 1.08
            ),
        },
        "series": [snv_series, indel_series],
    }


def shares_option(children: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Per child: paternal share of phased strict candidates, transitions."""
    rows = sorted(
        children,
        key=lambda c: (c.get("paternal_share") is None, c.get("paternal_share") or 0),
    )
    expected_pat = denovo.EXPECTED.get("paternal_share")
    expected_ti = denovo.EXPECTED.get("ti_fraction")

    def series(
        name: str, key: str, color: str, expected: Optional[float]
    ) -> Dict[str, Any]:
        s: Dict[str, Any] = {
            "name": name,
            "type": "scatter",
            "symbolSize": 9,
            "color": color,
            "data": [c.get(key) for c in rows],
        }
        if expected:
            s["markLine"] = {
                "silent": True,
                "symbol": "none",
                "lineStyle": {"color": color, "type": "dashed"},
                "label": {
                    "formatter": f"{name} (literature): {{c}}",
                    "position": "insideEndTop",
                },
                "data": [{"yAxis": expected}],
            }
        return s

    return {
        "title": {
            "text": "Paternal share and transitions",
            "subtext": "Paternal share of the phased strict candidates; transitions"
            " among strict SNVs. Dashed: literature",
            "left": "center",
            "textStyle": {"fontSize": 14},
        },
        "tooltip": {"trigger": "axis"},
        "legend": {"top": 62},
        "grid": {"top": 95, "left": 50, "right": 20, "bottom": 90},
        "xAxis": {
            "type": "category",
            "data": [_label(c) for c in rows],
            "axisLabel": {"rotate": 60, "fontSize": 10},
        },
        "yAxis": {"type": "value", "min": 0, "max": 1},
        "series": [
            series("Paternal share", "paternal_share", "#1565c0", expected_pat),
            series("Transitions", "ti_fraction", "#6a1b9a", expected_ti),
        ],
    }


def flags_option(flags: List[Dict[str, Any]], n_candidates: int) -> Dict[str, Any]:
    """Share of the candidates each flag marks, and of those it alone marks."""
    rows = list(reversed(flags))
    return {
        "title": {
            "text": f"Flags over {n_candidates:,} candidates",
            "subtext": "Light: the candidates the flag marks; dark: those it alone marks",
            "left": "center",
            "textStyle": {"fontSize": 14},
        },
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "shadow"}},
        "grid": {"top": 60, "left": 130, "right": 30, "bottom": 30},
        "xAxis": {"type": "value", "name": "candidates"},
        "yAxis": {"type": "category", "data": [r["label"] for r in rows]},
        "series": [
            {
                "name": "Marked",
                "type": "bar",
                "barGap": "-100%",
                "data": [
                    {
                        "value": r["n"],
                        "itemStyle": {
                            "color": r["FlagBadges"][0]["color"],
                            "opacity": 0.35,
                        },
                    }
                    for r in rows
                ],
            },
            {
                "name": "Only flag",
                "type": "bar",
                "data": [
                    {
                        "value": r["alone"],
                        "itemStyle": {"color": r["FlagBadges"][0]["color"]},
                    }
                    for r in rows
                ],
            },
        ],
    }
