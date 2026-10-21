"""Static SVG plots for the family report.

The report is one self-contained HTML document, shown in an iframe and printed
by the browser. ECharts and Plotly draw in JavaScript, after the page loads,
so a chart can be missing from the print. These builders return SVG strings
instead. They are pure, and every piece of text goes through ``html.escape``.
"""

import html
import math
import random
from dataclasses import dataclass
from statistics import NormalDist
from typing import Dict, List, Optional, Sequence, Tuple

from .ancestry import get_region_color

_NORMAL = NormalDist()
_FONT = "font-family:Helvetica,Arial,sans-serif"


def _e(text: object) -> str:
    return html.escape(str(text), quote=True)


def percentile(z: float) -> float:
    """The standard normal distribution's value at ``z``, in percent."""
    return 100.0 * _NORMAL.cdf(z)


def ordinal(pct: float) -> str:
    """``"69th"``; below 1 and above 99 as ``"<1st"`` / ``">99th"``."""
    if pct < 1:
        return "<1st"
    if pct > 99:
        return ">99th"
    n = int(round(pct))
    suffix = (
        "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    )
    return f"{n}{suffix}"


@dataclass(frozen=True)
class Marker:
    """One person on a plot."""

    label: str
    value: float
    color: str


# ---------------------------------------------------------------------------
# PGS: one small normal distribution per trait
# ---------------------------------------------------------------------------

_Z_MAX = 4.0


def _label_rows(xs: Sequence[float], min_gap: float) -> List[int]:
    """A label row per marker (by x order) so no two labels in a row overlap."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    last_x: List[float] = []
    rows = [0] * len(xs)
    for i in order:
        for row, x in enumerate(last_x):
            if xs[i] - x >= min_gap:
                last_x[row] = xs[i]
                rows[i] = row
                break
        else:
            rows[i] = len(last_x)
            last_x.append(xs[i])
    return rows


def pgs_distribution_svg(markers: Sequence[Marker], width: int = 295) -> str:
    """A standard normal curve with one marker per person at their z-score.

    The tails past two standard deviations are shaded. A z beyond the drawn
    range sits on the edge, with an arrow. Labels are stacked above the curve
    so close scores stay readable.
    """
    pad_l, pad_r = 10, 10
    label_h = 11
    rows = _label_rows(
        [max(-_Z_MAX, min(_Z_MAX, m.value)) for m in markers], min_gap=1.6
    )
    n_rows = (max(rows) + 1) if rows else 0
    top = 4 + n_rows * label_h
    curve_h = 52
    base = top + curve_h
    height = base + 26
    plot_w = width - pad_l - pad_r

    def sx(z: float) -> float:
        return pad_l + (z + _Z_MAX) / (2 * _Z_MAX) * plot_w

    peak = _NORMAL.pdf(0)

    def sy(z: float) -> float:
        return base - _NORMAL.pdf(z) / peak * curve_h

    steps = 120
    pts = [(sx(-_Z_MAX + 2 * _Z_MAX * i / steps), sy(-_Z_MAX + 2 * _Z_MAX * i / steps))
           for i in range(steps + 1)]  # fmt: skip
    curve = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)

    def tail(z0: float, z1: float) -> str:
        n = 30
        seg = [
            (sx(z0 + (z1 - z0) * i / n), sy(z0 + (z1 - z0) * i / n))
            for i in range(n + 1)
        ]
        body = " ".join(f"{x:.1f},{y:.1f}" for x, y in seg)
        return (
            f'<polygon points="{sx(z0):.1f},{base} {body} {sx(z1):.1f},{base}"'
            ' fill="#e5e7eb"/>'
        )

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"'
        f' viewBox="0 0 {width} {height}" style="{_FONT}">',
        tail(-_Z_MAX, -2),
        tail(2, _Z_MAX),
        f'<polyline points="{curve}" fill="none" stroke="#6b7280" stroke-width="1.2"/>',
        f'<line x1="{pad_l}" y1="{base}" x2="{width - pad_r}" y2="{base}"'
        ' stroke="#9ca3af" stroke-width="1"/>',
    ]
    for z in (-3, -2, -1, 0, 1, 2, 3):
        x = sx(z)
        parts.append(
            f'<line x1="{x:.1f}" y1="{base}" x2="{x:.1f}" y2="{base + 3}" stroke="#9ca3af"/>'
        )
        tick = f"{z:+d}" if z else "0"
        parts.append(
            f'<text x="{x:.1f}" y="{base + 12}" font-size="8" text-anchor="middle"'
            f' fill="#6b7280">{tick}</text>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{base + 21}" font-size="7" text-anchor="middle"'
            f' fill="#9ca3af">{_e(ordinal(percentile(z)))}</text>'
        )
    for m, row in zip(markers, rows):
        z = max(-_Z_MAX, min(_Z_MAX, m.value))
        x = sx(z)
        y = sy(z)
        label_y = 4 + (n_rows - row) * label_h - 2
        parts.append(
            f'<line x1="{x:.1f}" y1="{label_y + 2}" x2="{x:.1f}" y2="{base}"'
            f' stroke="{_e(m.color)}" stroke-width="1" stroke-dasharray="2,2"/>'
        )
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3.2" fill="{_e(m.color)}"'
            ' stroke="white" stroke-width="0.8"/>'
        )
        arrow = "◀ " if m.value < -_Z_MAX else ""
        arrow_r = " ▶" if m.value > _Z_MAX else ""
        anchor = (
            "start"
            if x < pad_l + 30
            else ("end" if x > width - pad_r - 30 else "middle")
        )
        parts.append(
            f'<text x="{x:.1f}" y="{label_y}" font-size="8" text-anchor="{anchor}"'
            f' fill="{_e(m.color)}">{_e(arrow + m.label + arrow_r)}</text>'
        )
    parts.append("</svg>")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Ancestry: PC1 x PC2 over the reference panel
# ---------------------------------------------------------------------------

Window = Tuple[Tuple[float, float], Tuple[float, float]]

#: Reference genomes drawn at most, sampled with a fixed seed so a report
#: printed twice draws the same points.
_MAX_REFERENCE_POINTS = 2500


def pca_svg(
    markers: Sequence[Tuple[str, float, float, str]],
    reference: Optional[Sequence[Tuple[float, float, str]]],
    window: Optional[Window],
    x_label: str = "PC1",
    y_label: str = "PC2",
    width: int = 290,
    height: int = 230,
) -> str:
    """The family's samples (``(label, x, y, color)``) over the reference panel.

    ``reference`` holds ``(x, y, region)`` points, drawn small and faded in
    their region's colour, with a legend of the regions shown. ``window``
    bounds the axes (see :func:`ancestry_plots.zoom_window`); without one the
    axes cover every point.
    """
    pad_l, pad_r, pad_t, pad_b = 34, 10, 10, 30
    ref = list(reference or [])
    if window is not None:
        (x0, x1), (y0, y1) = window
        ref = [p for p in ref if x0 <= p[0] <= x1 and y0 <= p[1] <= y1]
    if len(ref) > _MAX_REFERENCE_POINTS:
        ref = random.Random(0).sample(ref, _MAX_REFERENCE_POINTS)
    if window is None:
        xs = [m[1] for m in markers] + [p[0] for p in ref]
        ys = [m[2] for m in markers] + [p[1] for p in ref]
        if not xs:
            return ""
        span_x = (max(xs) - min(xs)) or 1e-3
        span_y = (max(ys) - min(ys)) or 1e-3
        window = (
            (min(xs) - 0.05 * span_x, max(xs) + 0.05 * span_x),
            (min(ys) - 0.05 * span_y, max(ys) + 0.05 * span_y),
        )
    (x0, x1), (y0, y1) = window
    x1 = x1 if x1 > x0 else x0 + 1e-3
    y1 = y1 if y1 > y0 else y0 + 1e-3
    plot_w = width - pad_l - pad_r
    plot_h = height - pad_t - pad_b

    def sx(v: float) -> float:
        return pad_l + (v - x0) / (x1 - x0) * plot_w

    def sy(v: float) -> float:
        return pad_t + (1 - (v - y0) / (y1 - y0)) * plot_h

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"'
        f' viewBox="0 0 {width} {height}" style="{_FONT}">',
        f'<rect x="{pad_l}" y="{pad_t}" width="{plot_w}" height="{plot_h}"'
        ' fill="none" stroke="#d1d5db"/>',
    ]
    regions: Dict[str, str] = {}
    for x, y, region in ref:
        color = get_region_color(region)
        regions.setdefault(region, color)
        parts.append(
            f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="1.3" fill="{_e(color)}"'
            ' fill-opacity="0.35"/>'
        )
    for label, x, y, color in markers:
        parts.append(
            f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="4" fill="{_e(color)}"'
            ' stroke="#111827" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{sx(x) + 6:.1f}" y="{sy(y) + 3:.1f}" font-size="9"'
            f' fill="#111827">{_e(label)}</text>'
        )
    parts.append(
        f'<text x="{pad_l + plot_w / 2:.1f}" y="{height - 16}" font-size="9"'
        f' text-anchor="middle" fill="#374151">{_e(x_label)}</text>'
    )
    parts.append(
        f'<text x="10" y="{pad_t + plot_h / 2:.1f}" font-size="9" text-anchor="middle"'
        f' fill="#374151" transform="rotate(-90 10 {pad_t + plot_h / 2:.1f})">'
        f"{_e(y_label)}</text>"
    )
    lx = pad_l
    for region, color in sorted(regions.items()):
        parts.append(
            f'<circle cx="{lx + 3}" cy="{height - 5}" r="3" fill="{_e(color)}"/>'
            f'<text x="{lx + 9}" y="{height - 2}" font-size="8" fill="#374151">'
            f"{_e(region)}</text>"
        )
        lx += 12 + 6 * len(region) + 8
    parts.append("</svg>")
    return "".join(parts)


def is_finite(value: object) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
    )
