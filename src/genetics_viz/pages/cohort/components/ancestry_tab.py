"""Ancestry tab component for the family page.

Shows the family's samples projected onto the versioned reference panel they
were scored against: PC1 x PC2 and PC3 x PC4, with the panel drawn as a
coloured background and a table of the predicted region per sample. A second
row repeats both plots zoomed to a window covering the family and every
reference sample of the regions predicted for it.
"""

import logging
import math
from typing import Any, Callable, Dict, List, Optional, Tuple

import polars as pl
from nicegui import ui

from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.ancestry import (
    NULL_VALUES,
    AncestryFiles,
    find_ancestry_files,
    get_ancestry_dir,
    get_reference_pcs_path,
    get_region_color,
    load_reference_pcs,
    probe_ancestry_data,
)
from genetics_viz.utils.tsv import read_tsv_or_none

logger = logging.getLogger(__name__)

__all__ = ["probe_ancestry_data", "render_ancestry_tab"]

#: Axis bounds for one plot: ``((x_min, x_max), (y_min, y_max))``.
Window = Tuple[Tuple[float, float], Tuple[float, float]]

#: The two principal component pairs plotted, in display order.
_PC_PAIRS = (("PC1", "PC2"), ("PC3", "PC4"))

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

_REGION_TABLE_COLUMNS = (
    ("IID", "Sample"),
    ("region", "Region"),
    ("population", "Population"),
    ("region_confidence", "Region conf."),
    ("population_confidence", "Pop. conf."),
    ("mean_knn_distance", "Mean kNN dist."),
    ("nearest_reference_distance", "Nearest ref. dist."),
    ("local_reference_density", "Local ref. density"),
    ("is_outlier", "Outlier"),
)


def _round(value: Any) -> Optional[float]:
    """Round a PC coordinate, returning None when it is not a number."""
    try:
        return round(float(value), _COORD_PRECISION)
    except (TypeError, ValueError):
        return None


def _axis_bounds(values: List[float]) -> Optional[Tuple[float, float]]:
    """Pad the extent of ``values`` into an axis range.

    The padded bounds are then rounded outwards to a step derived from the
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


def _zoom_window(
    rows: List[Dict],
    reference: Optional[pl.DataFrame],
    assigned_regions: List[str],
    x_pc: str,
    y_pc: str,
) -> Optional[Window]:
    """Compute the zoom window for one PC pair.

    The window covers every family sample plus every reference sample of the
    regions predicted for the family, so the family can be read against its own
    ancestry cluster instead of the whole panel. With no predicted regions - an
    ``ancestry.tsv`` that is missing or holds no region - it falls back to the
    family's own extent.

    Returns ``None`` when there is nothing to bound.
    """
    xs = [v for r in rows if (v := _round(r.get(x_pc))) is not None]
    ys = [v for r in rows if (v := _round(r.get(y_pc))) is not None]

    if reference is not None and assigned_regions:
        assigned = reference.filter(pl.col("region").is_in(assigned_regions))
        xs.extend(v for v in assigned[x_pc].to_list() if v is not None)
        ys.extend(v for v in assigned[y_pc].to_list() if v is not None)

    x_bounds, y_bounds = _axis_bounds(xs), _axis_bounds(ys)
    if x_bounds is None or y_bounds is None:
        return None
    return x_bounds, y_bounds


def _reference_series(
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
            x, y = _round(raw_x), _round(raw_y)
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


def _family_series(rows: List[Dict], x_pc: str, y_pc: str) -> Dict:
    """Build the foreground scatter series for the family's own samples."""
    data = []
    for row in rows:
        x, y = _round(row.get(x_pc)), _round(row.get(y_pc))
        if x is None or y is None:
            continue
        region = row.get("region") or ""
        data.append(
            {
                "name": row.get("IID", ""),
                "value": [
                    x,
                    y,
                    region or "—",
                    row.get("population") or "—",
                    row.get("region_confidence"),
                ],
                "itemStyle": {
                    "color": get_region_color(region),
                    "borderColor": "#d32f2f" if row.get("is_outlier") else "#1f2937",
                    "borderWidth": 2 if row.get("is_outlier") else 1.5,
                },
            }
        )
    return {
        "name": "Family",
        "type": "scatter",
        "data": data,
        "symbol": "diamond",
        "symbolSize": 14,
        "z": 10,
        "label": {
            "show": True,
            "position": "top",
            "fontSize": 10,
            "formatter": "{b}",
        },
        # Members of one family sit almost on top of each other in PC space,
        # so their labels would otherwise pile into an unreadable blob. The
        # tooltip still names every point.
        "labelLayout": {"hideOverlap": True},
    }


def _tooltip_formatter(x_pc: str, y_pc: str) -> str:
    """JS tooltip formatter showing the prediction for a family point.

    The reference background is ``silent``, so only the family series reaches
    this; the guard is there in case that ever changes.
    """
    return (
        "function (p) {"
        "  if (p.seriesName !== 'Family') { return ''; }"
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


def _plot_option(
    rows: List[Dict],
    reference: Optional[pl.DataFrame],
    x_pc: str,
    y_pc: str,
    window: Optional[Window] = None,
    title_suffix: str = "",
) -> Dict:
    """Build the ECharts option dict for one PC pair.

    ``window`` bounds both axes explicitly; without it the axes autoscale to
    the data. Explicit bounds matter for a zoomed plot - left to autoscale,
    ECharts would fit the filtered points and add its own ``scale`` padding, so
    the rendered window would not be the one that was computed.
    """
    series = (
        _reference_series(reference, x_pc, y_pc, window)
        if reference is not None
        else []
    )
    legend_names = [s["name"] for s in series]
    series.append(_family_series(rows, x_pc, y_pc))
    legend_names.append("Family")

    if window is None:
        x_axis_range: Dict[str, Any] = {"scale": True}
        y_axis_range: Dict[str, Any] = {"scale": True}
    else:
        (x0, x1), (y0, y1) = window
        x_axis_range = {"min": x0, "max": x1}
        y_axis_range = {"min": y0, "max": y1}

    return {
        "title": {
            "text": f"{x_pc} vs {y_pc}{title_suffix}",
            "left": "center",
            "textStyle": {"fontSize": 13, "fontWeight": "normal"},
        },
        # The ":" key prefix makes NiceGUI evaluate the value as a JS
        # function rather than passing it through as a string.
        "tooltip": {"trigger": "item", ":formatter": _tooltip_formatter(x_pc, y_pc)},
        "legend": {"data": legend_names, "bottom": 0, "type": "scroll"},
        # The grid leaves room below for the x-axis name *and* the legend
        # underneath it; a smaller bottom margin makes the two collide.
        "grid": {"top": 40, "bottom": 72, "left": 60, "right": 20},
        "xAxis": {
            "type": "value",
            "name": x_pc,
            "nameLocation": "middle",
            "nameGap": 25,
            **x_axis_range,
        },
        "yAxis": {
            "type": "value",
            "name": y_pc,
            "nameLocation": "middle",
            "nameGap": 45,
            **y_axis_range,
        },
        "series": series,
    }


def _render_region_table(rows: List[Dict]) -> None:
    """Render the predicted-region table for the selected members."""
    present = [col for col, _ in _REGION_TABLE_COLUMNS if any(col in r for r in rows)]
    if not present:
        return

    table_rows = []
    for row in rows:
        entry = {col: row.get(col) for col in present}
        if "is_outlier" in entry:
            entry["is_outlier"] = "yes" if row.get("is_outlier") else ""
        if "region" in entry:
            entry["region"] = row.get("region") or "—"
            entry["region_color"] = get_region_color(row.get("region"))
        table_rows.append(entry)

    headers = dict(_REGION_TABLE_COLUMNS)
    columns: List[Dict[str, Any]] = []
    for col in present:
        col_def: Dict[str, Any] = {
            "id": col,
            "header": headers[col],
            "sortable": True,
        }
        if col == "region":
            col_def["cellType"] = "badge"
            col_def["colorField"] = "region_color"
        elif col in (
            "region_confidence",
            "population_confidence",
            "mean_knn_distance",
            "nearest_reference_distance",
            "local_reference_density",
        ):
            col_def["cellType"] = "number"
            col_def["sorting"] = "numerical"
        columns.append(col_def)

    with ui.card().classes("w-full"):
        DataTable(
            columns=columns,
            rows=table_rows,
            row_key="IID",
            pagination={"rowsPerPage": 10},
        )


def _load_frames(files: AncestryFiles) -> Optional[List[Dict]]:
    """Read the PCs and join the ancestry predictions onto them.

    Returns row dicts, or ``None`` when the PCs hold no data. The join is a
    left join so the plots still render when ``ancestry.tsv`` is absent.
    """
    pcs = read_tsv_or_none(files.pcs_path, null_values=NULL_VALUES)
    if pcs is None or len(pcs) == 0:
        return None

    if files.ancestry_path is not None:
        try:
            predictions = read_tsv_or_none(files.ancestry_path, null_values=NULL_VALUES)
            if predictions is not None and len(predictions) > 0:
                overlap = (set(predictions.columns) & set(pcs.columns)) - {"IID"}
                if overlap:
                    predictions = predictions.drop(list(overlap))
                pcs = pcs.join(predictions, on="IID", how="left")
        except Exception:
            logger.warning("Failed to read %s", files.ancestry_path, exc_info=True)

    return pcs.to_dicts()


def render_ancestry_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the Ancestry tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    ancestry_dir = get_ancestry_dir(store.data_dir, family_id)
    if not ancestry_dir.exists():
        ui.label(f"No ancestry directory found at: {ancestry_dir}").classes(
            "text-gray-500 italic"
        )
        return

    files = find_ancestry_files(store.data_dir, family_id)
    if files is None or files.pcs_path is None:
        ui.label(f"No principal components found in: {ancestry_dir}").classes(
            "text-gray-500 italic"
        )
        return

    try:
        rows = _load_frames(files)
    except Exception as e:
        ui.label(f"Error reading ancestry files: {e}").classes("text-red-500 mt-4")
        logger.warning("Failed to read ancestry files for %s", family_id, exc_info=True)
        return

    if rows is None:
        ui.label("No principal components in file").classes("text-gray-500 italic")
        return

    reference = load_reference_pcs(files.bundle_version)

    # Provenance line, plus a note when several bundle versions are present.
    with ui.row().classes("items-center gap-2 mb-2"):
        ui.label(files.label).classes("text-xs text-gray-500")
        if files.other_variants:
            others = ", ".join(f"{b} / {f}" for b, f in files.other_variants)
            ui.icon("info_outline", color="grey").classes(
                "text-sm cursor-help"
            ).tooltip(f"Also present, not shown: {others}")

    if reference is None:
        reference_path = get_reference_pcs_path(files.bundle_version)
        detail = (
            "set 'ancestry_reference_dir' in the application config"
            if reference_path is None
            else f"expected at {reference_path}"
        )
        ui.label(
            f"Reference panel unavailable for bundle {files.bundle_version}"
            f" — {detail}. Showing the family's samples only."
        ).classes("text-orange-600 text-sm mb-2")

    @ui.refreshable
    def render_ancestry_content() -> None:
        selected = [r for r in rows if r.get("IID") in selected_members["value"]]
        if not selected:
            ui.label("No members selected").classes("text-gray-500 italic")
            return

        _render_region_table(selected)

        # Full-panel overview.
        with ui.row().classes("w-full gap-4 mt-2"):
            for x_pc, y_pc in _PC_PAIRS:
                ui.echart(_plot_option(selected, reference, x_pc, y_pc)).classes(
                    "flex-1 h-96"
                )

        # Zoomed view. Without a reference panel the overview plots already fit
        # the family alone, so a zoom row would just duplicate them.
        if reference is None:
            return

        assigned = sorted({r["region"] for r in selected if r.get("region")})
        if assigned:
            suffix = f" — zoom: {', '.join(assigned)}"
            caption = (
                "Zoomed to the family and every reference sample of "
                f"{', '.join(assigned)}."
            )
        else:
            suffix = " — zoom: family"
            caption = (
                "Zoomed to the family's own samples — no predicted region"
                " available to widen the window to."
            )

        ui.label(caption).classes("text-xs text-gray-500 mt-4")
        with ui.row().classes("w-full gap-4"):
            for x_pc, y_pc in _PC_PAIRS:
                window = _zoom_window(selected, reference, assigned, x_pc, y_pc)
                ui.echart(
                    _plot_option(
                        selected,
                        reference,
                        x_pc,
                        y_pc,
                        window=window,
                        title_suffix=suffix,
                    )
                ).classes("flex-1 h-96")

    render_ancestry_content()
    data_table_refreshers.append(render_ancestry_content.refresh)
