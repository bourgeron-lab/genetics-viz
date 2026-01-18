"""Waves validation page - displays samples from selected cohort pedigree."""

from collections import Counter
from typing import Dict, List

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.waves_loader import (
    get_wave_category,
    get_wave_color,
    load_waves_validations,
)
from genetics_viz.utils.data import get_data_store


@ui.page("/validation/waves")
def waves_validation_page() -> None:
    """Render the waves validation page showing samples from selected cohort."""
    create_header()

    try:
        store = get_data_store()
        validation_file = store.data_dir / "validations" / "waves.tsv"

        with ui.column().classes("w-full px-6 py-6"):
            ui.label("🌊 Waves Validation").classes(
                "text-3xl font-bold text-blue-900 mb-6"
            )

            # Cohort selector
            cohort_names = sorted(store.cohorts.keys())
            if not cohort_names:
                ui.label("No cohorts available").classes("text-gray-500 text-lg italic")
                return

            selected_cohort = {"value": cohort_names[0] if cohort_names else None}
# Filter state - all statuses selected by default
            all_statuses = ["TODO", "good", "low wave", "medium wave", "high wave"]
            filter_status: Dict[str, List[str]] = {"value": list(all_statuses)}

            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label("Select Cohort:").classes("text-lg font-semibold")
                cohort_select = (
                    ui.select(
                        options=cohort_names,
                        value=selected_cohort["value"],
                        on_change=lambda e: update_cohort(e.value),
                    )
                    .props("outlined dense")
                    .classes("w-64")
                )

            # Status filter button
            status_checkboxes: Dict[str, ui.checkbox] = {}

            with (
                ui.button(
                    "Filter by Status",
                    icon="filter_list",
                )
                .props("outline color=blue")
                .classes("mb-4")
            ):
                with ui.menu():
                    ui.label("Select Status:").classes("px-4 py-2 font-semibold text-sm")
                    ui.separator()

                    with ui.column().classes("p-2"):
                        # Quick select buttons
                        with ui.row().classes("gap-2 mb-2 flex-wrap"):

                            def select_all():
                                filter_status["value"] = list(all_statuses)
                                for cb in status_checkboxes.values():
                                    cb.value = True
                                refresh_content()

                            def select_none():
                  Prepare all table data first
                all_table_rows = []
                for sample_id in samples_list:
                    waves = waves_map.get(sample_id, [])
                    category = get_wave_category(waves)

                    # Format waves as colored badges
                    waves_display = (
                        ", ".join([str(w) for w in waves]) if waves else "None"
                    )

                    all_table_rows.append(
                        {
                            "sample_id": sample_id,
                            "validations": waves_display,
                            "count": len(waves),
                            "category": category,
                            "color": get_wave_color(category),
                        }
                    )

                # Apply filter
                filtered_table_rows = [
                    row
                    for row in all_table_rows
                    if row["category"] in filter_status["value"]
                ]

                # Calculate category counts (from filtered data)
                category_counts = Counter()
                for row in filtered_table_rows:
                    category_counts[row["category"]
                                "All",
                                on_click=select_all,
                            ).props("size=sm flat dense").classes("text-xs")
                            ui.button(
                                "None",
                                on_click=select_none,
                            ).props("size=sm flat dense").classes("text-xs")

                        ui.separator()

                        # Checkboxes for each status
                        for status in all_statuses:

                            def handle_change(st: str, val: bool) -> None:
                                if val:
                                    if st not in filter_status["value"]:
                                        filter_status["value"].append(st)
                                else:
                                    if st in filter_status["value"]:
                                        filter_status["value"].remove(st)
                                refresh_content()

                            color = get_wave_color(status)
                            with ui.row().classes("items-center gap-1"):
                                ui.badge("", color=color).classes("w-4 h-4")
                                status_checkboxes[status] = ui.checkbox(
                                    status,
                                    value=status in filter_status["value"],
                                    on_change=lambda e, s=status: handle_change(s, e.value),
                                ).classes("text-sm"    .classes("w-64")
                )

            # Containers for dynamic content
            pie_chart_container = ui.column().classes("w-full")
            table_container = ui.column().classes("w-full")

            def update_cohort(cohort_name: str):
                """Update th# Show filter count
                            if filter_status["value"] != all_statuses:
                                ui.label(
                                    f"Samples ({len(filtered_table_rows)} of {len(all_table_rows)})"
                                ).classes("text-xl font-semibold mb-4 text-blue-800")
                            else:
                                ui.label(f"Samples ({len(filtered_table_rows)})").classes(
                                    "text-xl font-semibold mb-4 text-blue-800"
                    refresh_content()

            def refresh_content():
                """Refresh the pie chart and table based on selected cohort."""
                if not selected_cohort["value"]:
                    return

                cohort = store.get_cohort(selected_cohort["value"])
                if not cohort:
                    return

                # Get all unique sample IDs from the pedigree
                sample_ids = set()
                for family in cohort.families.values():
                    for sample in family.samples:
                        sample_ids.add(sample.sample_id)

                samples_list = sorted(sample_ids)

                # Load validation data
                waves_map = load_waves_validations(validation_file)

                # Calculate category counts
                category_counts = Counter()
                for sample_id in samples_list:
                    waves = waves_map.get(sample_id)
                    category = get_wave_category(waves)
                    category_counts[category] += 1

                # Update pie chart
                pie_chart_container.clear()
                with pie_chart_container:
                    if samples_list:
                        with ui.card().classes("w-full mb-6 p-4"):
                            ui.label("Validation Status Overview").classes(
                                "text-xl font-semibold mb-4 text-blue-800"
                            )

                            # Prepare pie chart data
                            categories = [
                                "TODO",
                                "good",
                                "low wave",
                                "medium wave",
                                "high wave",
                            ]
                            pie_data = []
                            for category in categories:
                                count = category_counts.get(category, 0)
                                if count > 0:
                                    color = get_wave_color(category)
                                    pie_data.append(
                                        {
                                            "value": count,
                                            "name": f"{category} ({count})",
                                            "itemStyle": {"color": color},
                                        }
                                    )

                            echart_options = {
                                "tooltip": {"trigger": "item"},
                                "legend": {
                                    "orient": "vertical",
                                    "left": "left",
                                },
                                "series": [
                                    {
                                        "type": "pie",
                                        "radius": "50%",
                                        "data": pie_data,
                                        "emphasis": {
                                            "itemStyle": {
                                                "shadowBlur": 10,
                                                "shadowOffsetX": 0,
                                                "shadowColor": "rgba(0, 0, 0, 0.5)",
                                            }
                                        },
                                    }
                                ],
                            }
                            ui.echart(echart_options).classes("w-full h-80")

                # Update samples table
                table_container.clear()
                with table_container:
                    if not samples_list:
                        ui.label("No samples found in this cohort").classes(
                            "text-gray-500 text-lg italic"
                        )
                    else:
                        with ui.card().classes("w-full p-4"):
                            # Show filter count
                            if filter_status["value"] != all_statuses:
                                ui.label(
                                    f"Samples ({len(filtered_table_rows)} of {len(all_table_rows)})"
                                ).classes("text-xl font-semibold mb-4 text-blue-800")
                            else:
                                ui.label(f"Samples ({len(filtered_table_rows)})").classes(
                                    "text-xl font-semibold mb-4 text-blue-800"
                                )

                            # Prepare table data (already filtered)
                            table_rows = filtered_table_rows

                            # Table columns
                            columns = [
                                {
                                    "name": "sample_id",
                                    "label": "Sample ID",
                                    "field": "sample_id",
                                    "align": "left",
                                    "sortable": True,
                                },
                                {
                                    "name": "validations",
                                    "label": "Validations",
                                    "field": "validations",
                                    "align": "left",
                                },
                                {
                                    "name": "count",
                                    "label": "# Validations",
                                    "field": "count",
                                    "align": "center",
                                    "sortable": True,
                                },
                                {
                                    "name": "category",
                                    "label": "Status",
                                    "field": "category",
                                    "align": "center",
                                    "sortable": True,
                                },
                                {
                                    "name": "actions",
                                    "label": "Actions",
                                    "field": "sample_id",
                                    "align": "center",
                                },
                            ]

                            # Custom table slot for status with colors
                            table = ui.table(
                                columns=columns,
                                rows=table_rows,
                                row_key="sample_id",
                                pagination={"rowsPerPage": 20, "sortBy": "sample_id"},
                            ).classes("w-full")

                            # Add custom slot for category with color
                            table.add_slot(
                                "body-cell-category",
                                r"""
                                <q-td :props="props">
                                    <q-badge :color="props.row.color" :label="props.value" />
                                </q-td>
                                """,
                            )

                            # Add custom slot for actions (View on IGV button)
                            table.add_slot(
                                "body-cell-actions",
                                r"""
                                <q-td :props="props">
                                    <q-btn
                                        dense
                                        flat
                                        color="primary"
                                        icon="visibility"
                                        label="View on IGV"
                                        @click="$parent.$emit('view_igv', props.value)"
                                    />
                                </q-td>
                                """,
                            )

                            # Handle View on IGV button clicks
                            def handle_view_igv(e):
                                sample_id = e.args
                                ui.navigate.to(f"/validation/wave/{sample_id}")

                            table.on("view_igv", handle_view_igv)

            # Initial load
            refresh_content()

    except Exception as e:
        ui.label(f"Error loading waves validation page: {e}").classes(
            "text-red-600 text-lg"
        )
        import traceback

        ui.label(traceback.format_exc()).classes("text-xs text-gray-600 font-mono")
