"""Waves validation page - displays samples from selected cohort pedigree."""

import getpass
import json
from collections import Counter
from datetime import datetime
from typing import Dict, List

from nicegui import app as nicegui_app
from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.waves_loader import (
    get_wave_category,
    get_wave_color,
    get_wave_score_color,
    load_waves_validations,
    load_waves_validations_full,
    save_wave_validation,
)
from genetics_viz.utils.data import get_data_store


@ui.page("/validation/waves")
def waves_validation_page() -> None:
    """Render the waves validation page showing samples from selected cohort."""
    create_header()

    # Add IGV.js library at page level
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
    """)

    try:
        store = get_data_store()
        validation_file = store.data_dir / "validations" / "waves.tsv"

        # Serve data files for IGV.js
        nicegui_app.add_static_files("/data", str(store.data_dir))

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
                    ui.label("Select Status:").classes(
                        "px-4 py-2 font-semibold text-sm"
                    )
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
                                filter_status["value"] = []
                                for cb in status_checkboxes.values():
                                    cb.value = False
                                refresh_content()

                            ui.button(
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
                                    on_change=lambda e, s=status: handle_change(
                                        s, e.value
                                    ),
                                ).classes("text-sm")

            # Containers for dynamic content
            pie_chart_container = ui.column().classes("w-full")
            table_container = ui.column().classes("w-full")

            def update_cohort(cohort_name: str):
                """Update the display when cohort selection changes."""
                selected_cohort["value"] = cohort_name
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

                # Prepare all table data first
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
                    category_counts[row["category"]] += 1

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
                                ui.label(
                                    f"Samples ({len(filtered_table_rows)})"
                                ).classes("text-xl font-semibold mb-4 text-blue-800")

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
                                show_wave_dialog(sample_id)

                            table.on("view_igv", handle_view_igv)

            def show_wave_dialog(sample_id: str):
                """Show wave validation dialog for a sample."""
                # Check if bedgraph exists
                bedgraph_path = (
                    store.data_dir
                    / "samples"
                    / sample_id
                    / "sequences"
                    / f"{sample_id}.by1000.bedgraph.gz"
                )

                with (
                    ui.dialog().props("maximized") as dialog,
                    ui.card().classes("w-full h-full"),
                ):
                    with ui.column().classes("w-full h-full p-6"):
                        # Header with close button
                        with ui.row().classes(
                            "items-center justify-between w-full mb-4"
                        ):
                            ui.label(f"🌊 Wave Validation: {sample_id}").classes(
                                "text-2xl font-bold text-blue-900"
                            )
                            ui.button(
                                icon="close", on_click=lambda: dialog.close()
                            ).props("flat round")

                        if not bedgraph_path.exists():
                            with ui.card().classes("w-full p-6 bg-red-50"):
                                ui.label("⚠️ Bedgraph file not found").classes(
                                    "text-xl font-semibold text-red-800"
                                )
                                ui.label(f"Expected path: {bedgraph_path}").classes(
                                    "text-gray-600 text-sm font-mono"
                                )
                        else:
                            # IGV.js viewer
                            with ui.card().classes("w-full"):
                                igv_container = (
                                    ui.element("div")
                                    .classes("w-full")
                                    .style("height: 570px;")
                                )
                                igv_id = f"igv-{id(igv_container)}"
                                igv_container._props["id"] = igv_id

                                # Add IGV.js library if not already added
                                ui.add_head_html("""
                                    <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
                                """)

                                browser_var = f"igvBrowser_{igv_id.replace('-', '_')}"

                                # Build IGV config with bedgraph track
                                bedgraph_url = f"/data/samples/{sample_id}/sequences/{sample_id}.by1000.bedgraph.gz"
                                bedgraph_index_url = bedgraph_url + ".tbi"

                                igv_config = {
                                    "genome": "hg38",
                                    "locus": "chr1",
                                    "tracks": [
                                        {
                                            "name": f"{sample_id} Coverage",
                                            "type": "wig",
                                            "format": "bedgraph",
                                            "url": bedgraph_url,
                                            "indexURL": bedgraph_index_url,
                                            "height": 350,
                                            "color": "rgb(0, 0, 150)",
                                            "autoscale": False,
                                            "min": 0,
                                            "max": 80,
                                        }
                                    ],
                                }

                                # Store IGV config and ID for later initialization
                                igv_init_config = {
                                    "igv_id": igv_id,
                                    "browser_var": browser_var,
                                    "config": igv_config,
                                }

                            # Validation Form
                            validation_history_container = ui.column().classes(
                                "w-full mb-4"
                            )

                            def load_validation_history():
                                """Load and display validation history."""
                                waves_map_full = load_waves_validations_full(
                                    validation_file
                                )
                                validation_history_container.clear()

                                validations = waves_map_full.get(sample_id, [])

                                with validation_history_container:
                                    ui.label("Previous validations:").classes(
                                        "font-semibold mb-2"
                                    )
                                    if not validations:
                                        ui.label("No validations recorded yet").classes(
                                            "text-gray-500 text-sm italic"
                                        )
                                    else:
                                        with ui.card().classes("w-full p-2"):
                                            with ui.row().classes("gap-2 flex-wrap"):
                                                for validation in validations:
                                                    wave = validation["wave"]
                                                    user = validation["user"]
                                                    timestamp = validation["timestamp"]

                                                    wave_label = {
                                                        0: "good",
                                                        1: "low wave",
                                                        2: "medium wave",
                                                        3: "high wave",
                                                    }.get(wave, f"unknown ({wave})")

                                                    color = get_wave_score_color(wave)

                                                    # Format timestamp for display
                                                    try:
                                                        from datetime import datetime

                                                        dt = datetime.fromisoformat(
                                                            timestamp
                                                        )
                                                        formatted_time = dt.strftime(
                                                            "%Y-%m-%d %H:%M:%S"
                                                        )
                                                    except:
                                                        formatted_time = timestamp

                                                    tooltip_text = f"User: {user}\nDate: {formatted_time}"

                                                    ui.badge(
                                                        wave_label, color=color
                                                    ).classes(
                                                        "text-sm px-3 py-1"
                                                    ).tooltip(tooltip_text)

                            with ui.card().classes("w-full p-4 mb-4"):
                                with ui.column().classes("gap-4"):
                                    default_user = getpass.getuser()

                                    with ui.row().classes("items-center gap-4 w-full"):
                                        ui.label("User:").classes("font-semibold")
                                        user_input = (
                                            ui.input("Username")
                                            .props("outlined dense")
                                            .classes("w-48")
                                        )
                                        user_input.value = default_user

                                        ui.label("Wave Score:").classes(
                                            "font-semibold ml-4"
                                        )
                                        wave_select = (
                                            ui.select(
                                                options={
                                                    0: "0 - good",
                                                    1: "1 - low wave",
                                                    2: "2 - medium wave",
                                                    3: "3 - high wave",
                                                },
                                                value=0,
                                            )
                                            .props("outlined dense")
                                            .classes("w-56")
                                        )

                                        ui.button(
                                            "Save Validation",
                                            icon="save",
                                            on_click=lambda: save_validation_and_refresh(),
                                        ).props("color=blue").classes("ml-4")

                                    def save_validation_and_refresh():
                                        """Save a wave validation."""
                                        user = user_input.value.strip()
                                        wave = wave_select.value

                                        if not user:
                                            ui.notify(
                                                "Please enter a username",
                                                type="warning",
                                            )
                                            return

                                        if wave is None or not isinstance(wave, int):
                                            ui.notify(
                                                "Please select a valid wave score",
                                                type="warning",
                                            )
                                            return

                                        # Generate timestamp
                                        timestamp = datetime.now().isoformat()

                                        try:
                                            save_wave_validation(
                                                validation_file,
                                                sample_id,
                                                user,
                                                wave,
                                                timestamp,
                                            )

                                            ui.notify(
                                                f"Validation saved: {sample_id} - Wave {wave}",
                                                type="positive",
                                            )

                                            # Reload validation history
                                            load_validation_history()

                                            # Refresh the main table to show updated data
                                            refresh_content()

                                        except Exception as e:
                                            ui.notify(
                                                f"Error saving validation: {e}",
                                                type="negative",
                                            )
                                            import traceback

                                            print(traceback.format_exc())

                            load_validation_history()

                dialog.open()

                # Initialize IGV after dialog is open
                if bedgraph_path.exists():

                    def init_igv():
                        ui.run_javascript(
                            f"""
                            var igvDiv = document.getElementById("{igv_init_config["igv_id"]}");
                            console.log("Trying to initialize IGV, div found:", igvDiv !== null);
                            console.log("IGV library loaded:", typeof igv !== 'undefined');
                            if (igvDiv && typeof igv !== 'undefined') {{
                                igv.createBrowser(igvDiv, {json.dumps(igv_init_config["config"])})
                                    .then(function(browser) {{
                                        window.{igv_init_config["browser_var"]} = browser;
                                        console.log("IGV browser created in dialog");
                                    }})
                                    .catch(function(error) {{
                                        console.error("Error creating IGV browser:", error);
                                    }});
                            }} else {{
                                console.error("IGV container not found or igv not loaded");
                            }}
                        """
                        )

                    ui.timer(0.5, init_igv, once=True)

            # Initial load
            refresh_content()

    except Exception as e:
        ui.label(f"Error loading waves validation page: {e}").classes(
            "text-red-600 text-lg"
        )
        import traceback

        ui.label(traceback.format_exc()).classes("text-xs text-gray-600 font-mono")
