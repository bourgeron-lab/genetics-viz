"""Cohort-wide Ancestry tab for the cohort page.

Plots every member of the cohort on PC1 x PC2 over an optional reference-panel
background, with the distribution of predicted regions beside it. Everything
follows the filters applied to the individuals table.

Reads the cohort-level result files, ``<cohort>/ancestry/<cohort>.apgs_*.tsv``,
which carry an extra ``family_id`` column the per-family ones do not. Cohorts
whose pipeline run has not produced them get a disabled tab.
"""

import logging
from typing import Any, Callable, Dict, List, Optional, Tuple

from nicegui import ui

from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.ancestry import (
    NULL_VALUES,
    AncestryFiles,
    Variant,
    find_cohort_ancestry_files,
    get_cohort_ancestry_dir,
    get_region_color,
    get_reference_pcs_path,
    load_reference_pcs,
)
from genetics_viz.utils.ancestry_plots import (
    reference_series,
    round_coord,
    tooltip_formatter,
)
from genetics_viz.utils.tsv import read_tsv_or_none

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_ancestry_data", "render_cohort_ancestry_tab"]

#: The cohort plot shows one PC pair; the family tab is the place for PC3/PC4.
_X_PC, _Y_PC = "PC1", "PC2"

#: Name of the foreground series, and the label its legend entry carries.
_COHORT_SERIES = "Cohort"


def probe_cohort_ancestry_data(cohort: Any) -> bool:
    """Check whether cohort-level principal components exist for this cohort."""
    files = find_cohort_ancestry_files(cohort)
    return files is not None and files.pcs_path is not None


def _cohort_series(rows: List[Dict]) -> Dict:
    """Build the foreground scatter series for the cohort's own samples.

    Unlike the family tab this draws plain circles with no per-point labels:
    at 586 samples - let alone the 2676 of the largest cohort - labels are an
    unreadable smear. The tooltip still names every point.
    """
    data = []
    for row in rows:
        x, y = round_coord(row.get(_X_PC)), round_coord(row.get(_Y_PC))
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
                    "borderWidth": 1,
                },
            }
        )
    return {
        "name": _COHORT_SERIES,
        "type": "scatter",
        "data": data,
        "symbolSize": 7,
        "z": 10,
    }


def _plot_option(rows: List[Dict], ref_series: List[Dict]) -> Dict:
    """Build the ECharts option for the cohort PCA scatter.

    ``ref_series`` is passed in already built: it does not depend on the
    individuals filter, so the caller builds it once and hands over the same
    list on every refresh rather than re-walking ~6.7k reference rows per
    keystroke.
    """
    series = list(ref_series)
    legend_names = [s["name"] for s in series]
    series.append(_cohort_series(rows))
    legend_names.append(_COHORT_SERIES)

    return {
        "title": {
            "text": f"{_X_PC} vs {_Y_PC}",
            "left": "center",
            "textStyle": {"fontSize": 13, "fontWeight": "normal"},
        },
        # The ":" key prefix makes NiceGUI evaluate the value as a JS function.
        "tooltip": {
            "trigger": "item",
            ":formatter": tooltip_formatter(_X_PC, _Y_PC, _COHORT_SERIES),
        },
        "legend": {"data": legend_names, "bottom": 0, "type": "scroll"},
        # Room below for the x-axis name *and* the legend under it.
        "grid": {"top": 40, "bottom": 72, "left": 60, "right": 20},
        "xAxis": {
            "type": "value",
            "name": _X_PC,
            "nameLocation": "middle",
            "nameGap": 25,
            "scale": True,
        },
        "yAxis": {
            "type": "value",
            "name": _Y_PC,
            "nameLocation": "middle",
            "nameGap": 45,
            "scale": True,
        },
        "series": series,
    }


def _region_counts(rows: List[Dict]) -> List[Tuple[str, int]]:
    """Count samples per predicted region, most common first."""
    counts: Dict[str, int] = {}
    for row in rows:
        region = row.get("region") or "—"
        counts[region] = counts.get(region, 0) + 1
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))


def _render_region_summary(rows: List[Dict]) -> None:
    """Render the region distribution as a pie beside a count table."""
    counts = _region_counts(rows)
    if not counts:
        return
    total = sum(n for _, n in counts)

    with ui.row().classes("w-full gap-4 items-start"):
        pie_data = [
            {
                "value": n,
                "name": f"{region} ({n})",
                "itemStyle": {"color": get_region_color(region)},
            }
            for region, n in counts
        ]
        # No legend and no slice labels: the table beside the pie already
        # names every region with its colour, count and share, and at seven
        # slices the labels collide with each other and with the legend.
        ui.echart(
            {
                "tooltip": {"trigger": "item"},
                "series": [
                    {
                        "type": "pie",
                        "radius": ["35%", "75%"],
                        "data": pie_data,
                        "label": {"show": False},
                    }
                ],
            }
        ).classes("flex-1 h-64")

        table_rows = [
            {
                "region": region,
                "n": n,
                # Pre-formatted: the "number" cell type would render this as
                # "60.4000", which reads oddly for a percentage.
                "pct": f"{100 * n / total:.1f}%",
                "region_color": get_region_color(region),
            }
            for region, n in counts
        ]
        with ui.column().classes("flex-1").style("min-width: 0"):
            DataTable(
                columns=[
                    {
                        "id": "region",
                        "header": "Region",
                        "sortable": True,
                        "cellType": "badge",
                        "colorField": "region_color",
                    },
                    {
                        "id": "n",
                        "header": "N",
                        "sortable": True,
                        "sorting": "numerical",
                        "cellType": "number",
                    },
                    {
                        "id": "pct",
                        "header": "%",
                        "sortable": True,
                    },
                ],
                rows=table_rows,
                row_key="region",
                pagination={"rowsPerPage": 10},
            )


def _load_rows(files: AncestryFiles) -> Optional[List[Dict]]:
    """Read the cohort PCs with the ancestry predictions joined on.

    Left-joined, so the plot still renders when ``ancestry.tsv`` is absent -
    the points then fall back to the unknown-region colour.
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


def render_cohort_ancestry_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort-wide Ancestry tab.

    Args:
        store: DataStore instance
        cohort: Cohort instance
        filtered_state: Dict with 'individuals' key holding the filtered rows
        refreshers: List to append the refresh callback to
    """
    ancestry_dir = get_cohort_ancestry_dir(cohort)
    files = find_cohort_ancestry_files(cohort)
    if files is None or files.pcs_path is None:
        ui.label(f"No cohort-level ancestry results in: {ancestry_dir}").classes(
            "text-gray-500 italic"
        )
        return

    # Everything the filter touches is recomputed from these in-memory rows;
    # the TSVs are read once per selected variant, never per filter event.
    state: Dict[str, Any] = {
        "files": files,
        "rows": None,
        "ref_series": [],
        "reference_missing": False,
    }
    show_reference = {"value": True}

    def load_variant(variant: Optional[Variant] = None) -> None:
        """Read the frames for a variant and rebuild the reference series."""
        selected = find_cohort_ancestry_files(cohort, variant) or files
        state["files"] = selected
        try:
            state["rows"] = _load_rows(selected)
        except Exception:
            logger.warning(
                "Failed to read cohort ancestry for %s", cohort.name, exc_info=True
            )
            state["rows"] = None

        reference = load_reference_pcs(selected.bundle_version)
        state["reference_missing"] = reference is None
        # Independent of the individuals filter, so built once here and reused
        # on every refresh.
        state["ref_series"] = (
            reference_series(reference, _X_PC, _Y_PC) if reference is not None else []
        )

    load_variant()

    with ui.card().classes("w-full"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("Ancestry").classes("text-lg font-semibold text-blue-700")
            provenance = ui.label(state["files"].label).classes("text-xs text-gray-500")

        variants = state["files"].all_variants
        if len(variants) > 1:

            def on_variant(event: Any) -> None:
                bundle, tag = event.value.split(" / ", 1)
                load_variant((bundle, tag))
                provenance.text = state["files"].label
                render_content.refresh()

            ui.select(
                options=[f"{b} / {t}" for b, t in variants],
                value=f"{state['files'].bundle_version} / {state['files'].filter_tag}",
                label="Bundle / filters",
                on_change=on_variant,
            ).props("outlined dense").classes("w-full")

        ui.checkbox(
            "Show reference panel",
            value=show_reference["value"],
            on_change=lambda e: (
                show_reference.update({"value": e.value}),
                render_content.refresh(),
            ),
        )

        @ui.refreshable
        def render_content() -> None:
            rows = state["rows"]
            if rows is None:
                ui.label("No principal components in file").classes(
                    "text-gray-500 italic"
                )
                return

            if state["reference_missing"]:
                path = get_reference_pcs_path(state["files"].bundle_version)
                detail = (
                    "set 'ancestry_reference_dir' in the application config"
                    if path is None
                    else f"expected at {path}"
                )
                ui.label(
                    f"Reference panel unavailable for bundle"
                    f" {state['files'].bundle_version} — {detail}."
                ).classes("text-orange-600 text-sm")

            selected_ids = {
                ind["Sample ID"] for ind in filtered_state.get("individuals", [])
            }
            shown = [r for r in rows if r.get("IID") in selected_ids]
            if not shown:
                ui.label("No individuals match the current filters").classes(
                    "text-gray-500 italic"
                )
                return

            total = len(rows)
            caption = f"{len(shown)} individuals"
            if len(shown) < total:
                caption += f" (of {total} with ancestry results)"
            ui.label(caption).classes("text-xs text-gray-500")

            ref = state["ref_series"] if show_reference["value"] else []
            ui.echart(_plot_option(shown, ref)).classes("w-full h-96")
            _render_region_summary(shown)

        render_content()
        refreshers.append(render_content.refresh)
