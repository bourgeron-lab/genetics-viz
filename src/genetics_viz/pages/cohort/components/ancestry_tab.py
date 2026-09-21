"""Ancestry tab component for the family page.

Shows the family's samples projected onto the versioned reference panel they
were scored against: PC1 x PC2 and PC3 x PC4, with the panel drawn as a
coloured background and a table of the predicted region per sample. A second
row repeats both plots zoomed to a window covering the family and every
reference sample of the regions predicted for it.
"""

import logging
from typing import Any, Callable, Dict, List, Optional

import polars as pl
from nicegui import ui

from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.ancestry import (
    NULL_VALUES,
    AncestryFiles,
    Variant,
    find_ancestry_files,
    get_ancestry_dir,
    get_reference_pcs_path,
    get_region_color,
    load_reference_pcs,
    probe_ancestry_data,
)
from genetics_viz.utils.ancestry_plots import (
    Window,
    reference_series,
    round_coord,
    tooltip_formatter,
    zoom_window,
)
from genetics_viz.utils.tsv import read_tsv_or_none

logger = logging.getLogger(__name__)

__all__ = ["probe_ancestry_data", "render_ancestry_tab"]

#: The two principal component pairs plotted, in display order.
_PC_PAIRS = (("PC1", "PC2"), ("PC3", "PC4"))

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


def _family_series(rows: List[Dict], x_pc: str, y_pc: str) -> Dict:
    """Build the foreground scatter series for the family's own samples."""
    data = []
    for row in rows:
        x, y = round_coord(row.get(x_pc)), round_coord(row.get(y_pc))
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
        reference_series(reference, x_pc, y_pc, window) if reference is not None else []
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
        "tooltip": {"trigger": "item", ":formatter": tooltip_formatter(x_pc, y_pc)},
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

    # Reloaded when the variant selector changes; every render reads from here.
    state: Dict[str, Any] = {"files": files, "rows": None, "reference": None}

    def load_variant(variant: Optional[Variant] = None) -> bool:
        """Read the frames for a variant. False when they cannot be read."""
        selected = find_ancestry_files(store.data_dir, family_id, variant) or files
        state["files"] = selected
        try:
            state["rows"] = _load_frames(selected)
        except Exception:
            logger.warning(
                "Failed to read ancestry files for %s", family_id, exc_info=True
            )
            state["rows"] = None
            return False
        state["reference"] = load_reference_pcs(selected.bundle_version)
        return state["rows"] is not None

    if not load_variant():
        ui.label("No principal components in file").classes("text-gray-500 italic")
        return

    # Provenance line, with a selector when several variants are present.
    with ui.row().classes("items-center gap-2 mb-2"):
        provenance = ui.label(state["files"].label).classes("text-xs text-gray-500")

    if len(files.all_variants) > 1:

        def on_variant(event: Any) -> None:
            bundle, tag = event.value.split(" / ", 1)
            load_variant((bundle, tag))
            provenance.text = state["files"].label
            render_ancestry_content.refresh()

        ui.select(
            options=[f"{b} / {t}" for b, t in files.all_variants],
            value=f"{files.bundle_version} / {files.filter_tag}",
            label="Bundle / filters",
            on_change=on_variant,
        ).props("outlined dense").classes("w-64 mb-2")

    @ui.refreshable
    def render_ancestry_content() -> None:
        rows = state["rows"]
        reference = state["reference"]
        if rows is None:
            ui.label("No principal components in file").classes("text-gray-500 italic")
            return

        if reference is None:
            reference_path = get_reference_pcs_path(state["files"].bundle_version)
            detail = (
                "set 'ancestry_reference_dir' in the application config"
                if reference_path is None
                else f"expected at {reference_path}"
            )
            ui.label(
                f"Reference panel unavailable for bundle"
                f" {state['files'].bundle_version} — {detail}."
                f" Showing the family's samples only."
            ).classes("text-orange-600 text-sm mb-2")

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
                window = zoom_window(selected, reference, assigned, x_pc, y_pc)
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
