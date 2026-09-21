"""Shared ECharts builders for the ancestry PCA scatter plots.

Used by both the per-family Ancestry tab and the cohort-wide one. These are
pure functions returning option/series dicts - nothing here touches ``ui`` -
so they live in ``utils/`` rather than in one of the page components that
would otherwise have to import the other's privates.
"""

import math
from typing import Any, Dict, List, Optional, Tuple

import polars as pl

from .ancestry import get_region_color

#: Axis bounds for one plot: ``((x_min, x_max), (y_min, y_max))``.
Window = Tuple[Tuple[float, float], Tuple[float, float]]


#: Coordinates are rounded before being serialised - the reference panel
#: contributes ~6.7k points to each plot and dominates the payload.
_COORD_PRECISION = 5

#: Fraction of the span added on each side of a zoom window, so no point sits
#: exactly on an axis. Padding only ever grows the window, so the guarantee
#: that it covers the family and the assigned regions still holds.
_ZOOM_PADDING = 0.03

#: Pad applied when a zoom window's span is zero - a single sample, or samples
#: sharing a coordinate. Without it the axis would collapse to min == max.
_ZOOM_MIN_PAD = 1e-3

#: Axis bounds are rounded one digit finer than the point coordinates, so a
#: rounded point can never fall outside the window that contains it.
_BOUND_PRECISION = _COORD_PRECISION + 1


def round_coord(value: Any) -> Optional[float]:
    """Round a PC coordinate, returning None when it is not a number."""
    try:
        return round(float(value), _COORD_PRECISION)
    except (TypeError, ValueError):
        return None


def axis_bounds(values: List[float]) -> Optional[Tuple[float, float]]:
    """Pad the extent of ``values`` into an axis range.

    The padded bounds are rounded outwards to a step derived from the
    span. ECharts draws an explicit ``min``/``max`` as an end tick, so raw
    bounds would label the axis with something like ``-0.169274`` next to the
    round ticks beside it. Rounding outwards only grows the range, so every
    value stays inside it.
    """
    if not values:
        return None
    low, high = min(values), max(values)
    pad = (high - low) * _ZOOM_PADDING or _ZOOM_MIN_PAD
    low, high = low - pad, high + pad

    span = high - low
    if span <= 0:  # unreachable while pad > 0, but keeps log10 safe
        return round(low, _BOUND_PRECISION), round(high, _BOUND_PRECISION)
    step = 10 ** math.floor(math.log10(span)) / 5
    return (
        round(math.floor(low / step) * step, _BOUND_PRECISION),
        round(math.ceil(high / step) * step, _BOUND_PRECISION),
    )


def zoom_window(
    rows: List[Dict],
    reference: Optional[pl.DataFrame],
    assigned_regions: List[str],
    x_pc: str,
    y_pc: str,
) -> Optional[Window]:
    """Compute the zoom window for one PC pair.

    The window covers every sample in ``rows`` plus every reference sample of
    ``assigned_regions``, so the samples can be read against their own ancestry
    cluster instead of the whole panel. With no assigned regions - an
    ``ancestry.tsv`` that is missing or holds no region - it falls back to the
    extent of ``rows`` alone.

    Returns ``None`` when there is nothing to bound.
    """
    xs = [v for r in rows if (v := round_coord(r.get(x_pc))) is not None]
    ys = [v for r in rows if (v := round_coord(r.get(y_pc))) is not None]

    if reference is not None and assigned_regions:
        assigned = reference.filter(pl.col("region").is_in(assigned_regions))
        xs.extend(v for v in assigned[x_pc].to_list() if v is not None)
        ys.extend(v for v in assigned[y_pc].to_list() if v is not None)

    x_bounds, y_bounds = axis_bounds(xs), axis_bounds(ys)
    if x_bounds is None or y_bounds is None:
        return None
    return x_bounds, y_bounds


def reference_series(
    reference: pl.DataFrame,
    x_pc: str,
    y_pc: str,
    window: Optional[Window] = None,
) -> List[Dict]:
    """Build one scatter series per reference region.

    One series per region rather than one series overall, so the ECharts legend
    toggles regions for free. The background is ``silent`` so it does not steal
    tooltips from the family's own points.

    When a ``window`` is given, every region's points are kept but those
    outside the window are dropped: they would be clipped by the axes anyway,
    and dropping them here keeps the payload down. The regions that survive are
    the ones the legend offers, so a zoomed plot lists only what it can show.
    """
    if window is not None:
        (x0, x1), (y0, y1) = window
        reference = reference.filter(
            pl.col(x_pc).is_between(x0, x1) & pl.col(y_pc).is_between(y0, y1)
        )

    series: List[Dict] = []
    for (region,), group in reference.group_by(["region"], maintain_order=True):
        points = []
        for raw_x, raw_y in zip(group[x_pc].to_list(), group[y_pc].to_list()):
            x, y = round_coord(raw_x), round_coord(raw_y)
            if x is not None and y is not None:
                points.append([x, y])
        if not points:
            continue
        series.append(
            {
                "name": str(region),
                "type": "scatter",
                "data": points,
                "symbolSize": 4,
                "large": True,
                "largeThreshold": 2000,
                "silent": True,
                "itemStyle": {"color": get_region_color(str(region)), "opacity": 0.45},
            }
        )
    return series


def tooltip_formatter(x_pc: str, y_pc: str, series_name: str = "Family") -> str:
    """JS tooltip formatter showing the prediction for a foreground point.

    The reference background is ``silent``, so only ``series_name`` reaches
    this; the guard is there in case that ever changes.
    """
    return (
        "function (p) {"
        f"  if (p.seriesName !== '{series_name}') {{ return ''; }}"
        "  var v = p.value;"
        "  var conf = (v[4] === null || v[4] === undefined)"
        "    ? '\u2014' : Number(v[4]).toFixed(2);"
        "  return '<b>' + p.name + '</b><br/>'"
        "    + 'Region: ' + v[2] + ' (' + conf + ')<br/>'"
        "    + 'Population: ' + v[3] + '<br/>'"
        f"    + '{x_pc}: ' + v[0].toFixed(4) + '<br/>'"
        f"    + '{y_pc}: ' + v[1].toFixed(4);"
        "}"
    )
