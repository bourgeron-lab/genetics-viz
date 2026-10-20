"""Plot builders for the Sample QC tabs.

Pure functions returning ECharts option dicts - nothing here touches ``ui`` -
used by both the per-family tab and the cohort one. Every line drawn is a
threshold the pipeline used (``config/sample_qc.yaml``, ``roh.yaml``'s KING
cut-offs); the classes and flags themselves come from the pipeline.
"""

import html
import math
import random
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .pedigree_labels import SEX_COLORS
from .roh import KINSHIP_DEGREES
from .sample_qc import (
    HUB_MINIMUM,
    KARYOTYPE_POINTS,
    THRESHOLDS,
    CrossFamilyGraph,
    PairRow,
    relationship_info,
)

_PRECISION = 5
_BACKGROUND_COLOR = "rgba(120, 120, 120, 0.30)"
_MISMATCH_COLOR = "#c62828"
_TITLE_STYLE = {"fontSize": 12, "fontWeight": "normal"}

#: Pedigree sex as the samples table writes it -> SEX_COLORS code.
_SEX_CODE = {"male": "1", "female": "2"}

#: Families in the cross-family graph cycle through these.
_FAMILY_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]  # fmt: skip


def _round(value: Optional[float]) -> Optional[float]:
    return None if value is None else round(value, _PRECISION)


def _esc(text: Any) -> str:
    return html.escape(str(text), quote=True)


def sex_color(pedigree_sex: Optional[str]) -> str:
    return SEX_COLORS[_SEX_CODE.get(pedigree_sex or "", "-9")]


def _title(text: str) -> Dict[str, Any]:
    return {"text": text, "left": "center", "textStyle": _TITLE_STYLE}


def _line(axis: str, value: float, label: str, color: str = "#455a64") -> Dict:
    return {
        axis + "Axis": value,
        "name": label,
        "lineStyle": {"color": color, "type": "dashed", "width": 1},
        "label": {
            "formatter": label,
            "color": color,
            "fontSize": 10,
            "position": "insideEndTop" if axis == "y" else "insideStartTop",
        },
    }


def _tooltip_item() -> Dict[str, Any]:
    """Item tooltip reading the ``tip`` HTML each data item carries."""
    return {
        "trigger": "item",
        "confine": True,
        ":formatter": "p => (p.data && p.data.tip) || p.name || ''",
    }


def king_marklines() -> List[Dict]:
    """Dashed horizontal lines at the KING kinship cut-offs, labelled."""
    lines = []
    for entry in KINSHIP_DEGREES:
        if entry.get("degree") is None:
            continue
        lines.append(
            _line("y", float(entry["min"]), f"{entry['label']} ≥ {entry['min']:g}")
        )
    return lines


# ---------------------------------------------------------------------------
# Relatedness
# ---------------------------------------------------------------------------


def _pair_tip(p: PairRow, show_detail: bool) -> str:
    lines = [f"<b>{_esc(p.sample_a)} – {_esc(p.sample_b)}</b>"]
    if p.cross_family:
        lines.append(f"families {_esc(p.family_a)} / {_esc(p.family_b)}")
    if show_detail or p.declared:
        lines.append(f"pedigree: {_esc(relationship_info(p.expected)['label'])}")
    if show_detail:
        lines.append(f"observed: {_esc(relationship_info(p.observed)['label'])}")
    if p.kinship is not None:
        lines.append(f"KING kinship {p.kinship:.4f}")
    if p.ratio is not None:
        lines.append(f"IBS0/IBS2 {p.ratio:.4f}")
    if p.n_sites is not None:
        lines.append(f"{int(p.n_sites)} sites")
    if p.status == "mismatch":
        lines.append(f"<span style='color:{_MISMATCH_COLOR}'>mismatch</span>")
    return "<br>".join(lines)


def relatedness_scatter_option(
    pairs: Sequence[PairRow],
    *,
    title: str,
    show_detail: bool,
    label_mismatches: bool = True,
) -> Dict:
    """KING kinship against IBS0/IBS2: the somalier/peddy relatedness plane.

    Horizontal lines are the KING degree cut-offs, the vertical one the
    parent-offspring / full-sib split (IBS0/IBS2 below 0.01). Points are
    coloured by the relationship the pedigree *declares*, so an error is a
    point sitting in another class's band; mismatches are red-edged diamonds.
    Viewers get one colour: the declared classes of undeclared pairs would
    say nothing, and the parents' pair is a consanguinity reading.
    """
    by_series: Dict[str, List[Dict]] = {}
    mismatches: List[Dict] = []
    xs: List[float] = []
    ys: List[float] = []
    for p in pairs:
        if p.kinship is None or p.ratio is None:
            continue
        xs.append(p.ratio)
        ys.append(p.kinship)
        item = {
            "name": f"{p.sample_a} – {p.sample_b}",
            "value": [_round(p.ratio), _round(p.kinship)],
            "tip": _pair_tip(p, show_detail),
        }
        if p.status == "mismatch":
            mismatches.append(item)
            continue
        key = relationship_info(p.expected)["label"] if show_detail else "Pairs"
        by_series.setdefault(key, []).append(item)

    series: List[Dict] = []
    colors = {
        relationship_info(name)["label"]: relationship_info(name)["color"]
        for name in ("parent-offspring", "full-sibs", "first-degree", "half-sibs",
                     "second-degree", "third-degree", "unrelated", "duplicate")
    }  # fmt: skip
    marks = {
        "symbol": "none",
        "silent": True,
        "data": king_marklines()
        + [_line("x", THRESHOLDS.get("po_ibs0_ibs2_max", 0.01), "PO | FS")],
    }
    for i, (name, data) in enumerate(sorted(by_series.items())):
        series.append(
            {
                "name": name,
                "type": "scatter",
                "data": data,
                "symbolSize": 8,
                "itemStyle": {"color": colors.get(name, "#607d8b"), "opacity": 0.8},
                **({"markLine": marks} if i == 0 else {}),
            }
        )
    series.append(
        {
            "name": "Mismatch",
            "type": "scatter",
            "data": mismatches,
            "symbol": "diamond",
            "symbolSize": 13,
            "z": 10,
            "itemStyle": {
                "color": _MISMATCH_COLOR,
                "borderColor": "#1f2937",
                "borderWidth": 1,
            },
            "label": {
                "show": label_mismatches,
                "position": "right",
                "fontSize": 9,
                "formatter": "{b}",
            },
            "labelLayout": {"hideOverlap": True},
            **({"markLine": marks} if not by_series else {}),
        }
    )
    x_max = max(0.2, math.ceil(max(xs, default=0) * 1.1 * 100) / 100)
    y_min = min(-0.05, math.floor(min(ys, default=0) * 10) / 10)
    return {
        "title": _title(title),
        "tooltip": _tooltip_item(),
        "legend": {"bottom": 0, "type": "scroll", "textStyle": {"fontSize": 10}},
        "grid": {"top": 40, "bottom": 70, "left": 55, "right": 25},
        "xAxis": {
            "type": "value",
            "name": "IBS0 / IBS2",
            "nameLocation": "middle",
            "nameGap": 25,
            "min": 0,
            "max": x_max,
        },
        "yAxis": {
            "type": "value",
            "name": "KING kinship",
            "nameLocation": "middle",
            "nameGap": 40,
            "min": y_min,
            "max": 0.55,
        },
        "series": series,
    }


def kinship_matrix_option(
    members: Sequence[str], pairs: Sequence[PairRow], *, show_detail: bool
) -> Dict:
    """Members × members: observed KING kinship above the diagonal, expected below.

    One colour scale for both halves, so a pair reads as two cells of a
    different shade when observation and pedigree disagree; a mismatch is
    also outlined in red. Viewers' pairs are pre-filtered to the declared ones.
    """
    index = {m: i for i, m in enumerate(members)}
    data: List[Dict] = []
    for p in pairs:
        if p.sample_a not in index or p.sample_b not in index:
            continue
        i, j = sorted((index[p.sample_a], index[p.sample_b]))
        border = (
            {"borderColor": _MISMATCH_COLOR, "borderWidth": 3}
            if p.status == "mismatch"
            else {"borderColor": "#ffffff", "borderWidth": 1}
        )
        tip = _pair_tip(p, show_detail)
        if p.kinship is not None:
            # Upper triangle: column j, row i.
            data.append(
                {
                    "value": [j, i, _round(max(p.kinship, 0.0))],
                    "tip": "observed<br>" + tip,
                    "itemStyle": border,
                    "label": {"formatter": f"{p.kinship:.2f}"},
                }
            )
        if p.expected_kinship is not None:
            data.append(
                {
                    "value": [i, j, _round(p.expected_kinship)],
                    "tip": "pedigree<br>" + tip,
                    "itemStyle": border,
                    "label": {"formatter": f"{p.expected_kinship:.2f}"},
                }
            )
    return {
        "tooltip": _tooltip_item(),
        "grid": {"top": 45, "bottom": 80, "left": 80, "right": 20},
        "title": _title("Kinship: observed ↗ / pedigree ↙"),
        "xAxis": {
            "type": "category",
            "data": list(members),
            "position": "bottom",
            "axisLabel": {"rotate": 45, "fontSize": 10},
            "splitArea": {"show": True},
        },
        "yAxis": {
            "type": "category",
            "data": list(members),
            "inverse": True,
            "axisLabel": {"fontSize": 10},
            "splitArea": {"show": True},
        },
        "visualMap": {
            "min": 0,
            "max": 0.5,
            "calculable": False,
            "orient": "horizontal",
            "left": "center",
            "bottom": 0,
            "itemHeight": 120,
            "text": ["0.5", "0"],
            "inRange": {"color": ["#f5f5f5", "#90caf9", "#1565c0", "#0d1b4c"]},
        },
        "series": [
            {
                "type": "heatmap",
                "data": data,
                "label": {"show": True, "fontSize": 10},
                "emphasis": {"itemStyle": {"shadowBlur": 6}},
            }
        ],
    }


# ---------------------------------------------------------------------------
# Sex and karyotype
# ---------------------------------------------------------------------------


def _sample_tip(p: Dict[str, Any], show_detail: bool) -> str:
    lines = [f"<b>{_esc(p['name'])}</b>"]
    if p.get("family"):
        lines.append(f"family {_esc(p['family'])}")
    lines.append(f"pedigree sex: {_esc(p.get('pedigree_sex') or 'unknown')}")
    if show_detail:
        lines.append(f"inferred: {_esc(p.get('inferred_sex') or 'unknown')}")
        lines.append(f"karyotype: {_esc(p.get('karyotype') or 'unknown')}")
    for label, key in (("X copies", "x"), ("Y copies", "y"), ("X het", "het")):
        if p.get(key) is not None:
            lines.append(f"{label} {p[key]:.3f}")
    if p.get("evidence"):
        lines.append(f"from {_esc(p['evidence'])}")
    return "<br>".join(lines)


def _member_series(
    points: Sequence[Dict[str, Any]],
    x_key: str,
    y_key: str,
    show_detail: bool,
) -> List[Dict]:
    """The highlighted samples, coloured by pedigree sex; flagged ones as diamonds."""
    normal, flagged = [], []
    for p in points:
        x, y = p.get(x_key), p.get(y_key)
        if x is None or y is None:
            continue
        item = {
            "name": p["name"],
            "value": [_round(x), _round(y)],
            "tip": _sample_tip(p, show_detail),
            "itemStyle": {
                "color": sex_color(p.get("pedigree_sex")),
                "borderColor": "#1f2937",
                "borderWidth": 1,
            },
        }
        (flagged if p.get("flagged") else normal).append(item)
    label = {
        "show": True,
        "position": "right",
        "fontSize": 9,
        "formatter": "{b}",
    }
    return [
        {
            "name": "Samples",
            "type": "scatter",
            "data": normal,
            "symbolSize": 10,
            "z": 5,
            "label": {**label, "show": len(normal) <= 12},
            "labelLayout": {"hideOverlap": True},
        },
        {
            "name": "Flagged",
            "type": "scatter",
            "data": flagged,
            "symbol": "diamond",
            "symbolSize": 15,
            "z": 10,
            "itemStyle": {"borderColor": _MISMATCH_COLOR, "borderWidth": 2},
            "label": label,
            "labelLayout": {"hideOverlap": True},
        },
    ]


def _background_series(points: Iterable[Tuple[float, float]], **extra: Any) -> Dict:
    return {
        "name": "Cohort",
        "type": "scatter",
        "data": [[_round(x), _round(y)] for x, y in points],
        "symbolSize": 5,
        "itemStyle": {"color": _BACKGROUND_COLOR},
        "large": True,
        "silent": True,
        **extra,
    }


def sex_scatter_option(
    points: Sequence[Dict[str, Any]],
    *,
    background: Sequence[Tuple[float, float]] = (),
    title: str = "Sex chromosome copy numbers",
    show_detail: bool,
) -> Dict:
    """X copies against Y copies, over shaded boxes around each karyotype.

    The boxes are the pipeline's tolerance around a whole copy number: a
    sample outside every box is ``ambiguous`` (mosaic or artefact). Points are
    coloured by *pedigree* sex, so a blue point in the XX box is a mismatch
    at a glance.

    Args:
        points: ``{"name", "x", "y", "pedigree_sex", "flagged", ...}`` per sample.
        background: ``(x, y)`` of the rest of the cohort, in grey.
    """
    tol = THRESHOLDS.get("sex_tolerance", 0.25)
    areas = [
        [
            {
                "name": name,
                "xAxis": x - tol,
                "yAxis": y - tol,
                "itemStyle": {
                    "color": (
                        "rgba(46,125,50,0.10)"
                        if name in ("XX", "XY")
                        else "rgba(106,27,154,0.10)"
                    )
                },
                "label": {
                    "color": "#455a64",
                    "fontSize": 10,
                    "position": "insideTopLeft",
                },
            },
            {"xAxis": x + tol, "yAxis": y + tol},
        ]
        for name, (x, y) in KARYOTYPE_POINTS.items()
    ]
    xs = [p["x"] for p in points if p.get("x") is not None] + [x for x, _ in background]
    ys = [p["y"] for p in points if p.get("y") is not None] + [y for _, y in background]
    x_max = max(3.4, math.ceil(max(xs, default=0) + 0.3))
    y_max = max(2.4, math.ceil(max(ys, default=0) + 0.3))
    series = [
        _background_series(background, markArea={"silent": True, "data": areas}),
        *_member_series(points, "x", "y", show_detail),
    ]
    return {
        "title": _title(title),
        "tooltip": _tooltip_item(),
        "grid": {"top": 40, "bottom": 45, "left": 50, "right": 20},
        "xAxis": {
            "type": "value",
            "name": "X copies",
            "nameLocation": "middle",
            "nameGap": 25,
            "min": 0,
            "max": x_max,
        },
        "yAxis": {
            "type": "value",
            "name": "Y copies",
            "nameLocation": "middle",
            "nameGap": 30,
            "min": 0,
            "max": y_max,
        },
        "series": series,
    }


def x_het_option(
    points: Sequence[Dict[str, Any]],
    *,
    background: Sequence[Tuple[float, float]] = (),
    show_detail: bool,
) -> Dict:
    """chrX heterozygosity against X copies.

    One X should be almost never heterozygous (below the upper line), two or
    more often (above the lower one). A sample in the wrong band is the
    pipeline's ``X_HET_INCONSISTENT``: contamination, or an X copy number
    the coverage got wrong. An XXY sample sits with the females here.
    """
    one_max = THRESHOLDS.get("x_het_one_x_max", 0.25)
    two_min = THRESHOLDS.get("x_het_two_x_min", 0.1)
    areas = [
        [
            {
                "name": "expected for 1 X",
                "xAxis": 0,
                "yAxis": 0,
                "itemStyle": {"color": "rgba(30,136,229,0.08)"},
                "label": {
                    "fontSize": 10,
                    "color": "#455a64",
                    "position": "insideTopLeft",
                },
            },
            {"xAxis": 1.5, "yAxis": one_max},
        ],
        [
            {
                "name": "expected for 2+ X",
                "xAxis": 1.5,
                "yAxis": two_min,
                "itemStyle": {"color": "rgba(216,27,96,0.08)"},
                "label": {
                    "fontSize": 10,
                    "color": "#455a64",
                    "position": "insideTopRight",
                },
            },
            {"xAxis": "max", "yAxis": 1},
        ],
    ]
    xs = [p["x"] for p in points if p.get("x") is not None] + [x for x, _ in background]
    return {
        "title": _title("chrX heterozygosity"),
        "tooltip": _tooltip_item(),
        "grid": {"top": 40, "bottom": 45, "left": 50, "right": 20},
        "xAxis": {
            "type": "value",
            "name": "X copies",
            "nameLocation": "middle",
            "nameGap": 25,
            "min": 0,
            "max": max(3.4, math.ceil(max(xs, default=0) + 0.3)),
        },
        "yAxis": {
            "type": "value",
            "name": "X het / (het + hom alt)",
            "nameLocation": "middle",
            "nameGap": 35,
            "min": 0,
            "max": 1,
        },
        "series": [
            _background_series(background, markArea={"silent": True, "data": areas}),
            *_member_series(points, "x", "het", show_detail),
        ],
    }


# ---------------------------------------------------------------------------
# Sample quality
# ---------------------------------------------------------------------------


def metric_strip_option(
    points: Sequence[Dict[str, Any]],
    *,
    background: Sequence[float] = (),
    title: str,
    axis_name: str,
    lines: Sequence[Tuple[float, str, str]] = (),
) -> Dict:
    """One metric as a horizontal strip: the cohort jittered, samples as diamonds.

    Args:
        points: ``{"name", "value", "flagged"}`` per highlighted sample.
        background: Values of the rest of the cohort, in grey.
        lines: ``(value, label, colour)`` vertical lines (thresholds, fences).
    """
    rng = random.Random(0)
    bg = [[_round(v), round(rng.uniform(-1, 1), 3)] for v in background]
    n = len(points)
    data = []
    for i, p in enumerate(points):
        if p.get("value") is None:
            continue
        y = 0.0 if n == 1 else -0.7 + 1.4 * i / (n - 1)
        data.append(
            {
                "name": p["name"],
                "value": [_round(p["value"]), round(y, 3)],
                "tip": f"<b>{_esc(p['name'])}</b><br>{_esc(axis_name)}: {p['value']:.4f}",
                "itemStyle": {
                    "color": _MISMATCH_COLOR if p.get("flagged") else "#1565c0",
                    "borderColor": "#1f2937",
                    "borderWidth": 1,
                },
            }
        )
    values = [v for v in background] + [
        p["value"] for p in points if p.get("value") is not None
    ]
    top = max(values + [v for v, _, _ in lines], default=0.0)
    return {
        "title": _title(title),
        "tooltip": _tooltip_item(),
        "grid": {"top": 40, "bottom": 40, "left": 20, "right": 30},
        "xAxis": {
            "type": "value",
            "name": axis_name,
            "nameLocation": "middle",
            "nameGap": 25,
            "min": 0,
            "max": _nice_max(top),
            "axisLine": {"onZero": False},
        },
        "yAxis": {"type": "value", "show": False, "min": -1.3, "max": 1.3},
        "series": [
            _background_series(
                [(x, y) for x, y in bg],
                markLine={
                    "symbol": "none",
                    "silent": True,
                    "data": [_line("x", v, label, color) for v, label, color in lines],
                },
            ),
            {
                "name": "Samples",
                "type": "scatter",
                "data": data,
                "symbol": "diamond",
                "symbolSize": 13,
                "z": 10,
                "label": {
                    "show": n <= 12,
                    "position": "top",
                    "fontSize": 9,
                    "formatter": "{b}",
                },
                "labelLayout": {"hideOverlap": True},
            },
        ],
    }


def _nice_max(top: float) -> float:
    """An axis end a little above ``top``, on two significant digits."""
    if top <= 0:
        return 1.0
    magnitude = 10 ** math.floor(math.log10(top * 1.15))
    return math.ceil(top * 1.15 / magnitude) * magnitude


# ---------------------------------------------------------------------------
# Across families
# ---------------------------------------------------------------------------


def family_color(family_id: str, order: Dict[str, int]) -> str:
    return _FAMILY_PALETTE[
        order.setdefault(family_id, len(order)) % len(_FAMILY_PALETTE)
    ]


def cross_family_graph_option(
    graph: CrossFamilyGraph,
    *,
    show_detail: bool,
    highlight: Optional[Set[Tuple[str, str]]] = None,
    hub_minimum: int = HUB_MINIMUM,
) -> Dict:
    """The cross-family related pairs as a force-directed graph of samples.

    Node size grows with the number of relatives outside the family, so a
    sample related to many families - contamination or population structure
    far more often than real relatives - stands out as a hub. Edges are
    coloured by the observed relationship (curators) and nodes by family.
    """
    order: Dict[str, int] = {}
    nodes = []
    ids: Dict[Tuple[str, str], str] = {}
    for key, node in sorted(graph.nodes.items(), key=lambda kv: -kv[1].degree):
        node_id = f"{node.family_id}|{node.sample}"
        ids[key] = node_id
        is_hub = node.degree >= hub_minimum
        shown = highlight is None or key in highlight
        nodes.append(
            {
                "id": node_id,
                "name": node.sample,
                "value": node.degree,
                "symbolSize": 8 + 6 * math.sqrt(node.degree - 1),
                "itemStyle": {
                    "color": family_color(node.family_id, order),
                    "opacity": 1.0 if shown else 0.35,
                    "borderColor": _MISMATCH_COLOR if is_hub else "#ffffff",
                    "borderWidth": 2 if is_hub else 1,
                },
                "label": {"show": is_hub},
                "tip": (
                    f"<b>{_esc(node.sample)}</b><br>family {_esc(node.family_id)}"
                    f"<br>{node.degree} relative(s) in other families"
                    + (
                        "<br><i>hub: check contamination / population</i>"
                        if is_hub
                        else ""
                    )
                ),
            }
        )
    links = []
    for p in graph.edges:
        info = relationship_info(p.observed)
        tip = f"<b>{_esc(p.sample_a)} – {_esc(p.sample_b)}</b>"
        if show_detail:
            tip += f"<br>{_esc(info['label'])}"
            if p.kinship is not None:
                tip += f", KING {p.kinship:.4f}"
        links.append(
            {
                "source": ids[(p.family_a, p.sample_a)],
                "target": ids[(p.family_b, p.sample_b)],
                "tip": tip,
                "lineStyle": {
                    "color": info["color"] if show_detail else "#9e9e9e",
                    "width": 1 + 8 * max(p.kinship or 0.0, 0.0),
                    "opacity": 0.8,
                },
            }
        )
    return {
        "tooltip": _tooltip_item(),
        "series": [
            {
                "type": "graph",
                "layout": "force",
                "roam": True,
                "draggable": True,
                "data": nodes,
                "links": links,
                "force": {"repulsion": 80, "edgeLength": [20, 70], "gravity": 0.08},
                "label": {"position": "right", "fontSize": 10, "formatter": "{b}"},
                "labelLayout": {"hideOverlap": True},
                "emphasis": {"focus": "adjacency"},
            }
        ],
    }
