"""Validation all page - displays all validations from validations/snvs.tsv."""

import csv
from typing import Any, Dict, List

from nicegui import app as nicegui_app
from nicegui import ui

from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.header import create_header
from genetics_viz.components.tables import VALIDATION_TABLE_SLOT
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.data import get_data_store


@ui.page("/validation/all")
def validation_all_page() -> None:
    """Render all validations from validations/snvs.tsv."""
    create_header()

    # Add IGV.js library at page level
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
    """)

    try:
        store = get_data_store()
        validation_file = store.data_dir / "validations" / "snvs.tsv"

        # Serve data files for IGV.js
        nicegui_app.add_static_files("/data", str(store.data_dir))

        with ui.column().classes("w-full px-6 py-6"):
            # Title
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label("📋 All Validations").classes(
                    "text-3xl font-bold text-blue-900"
                )

            if not validation_file.exists():
                ui.label("No validations found").classes("text-gray-500 text-lg italic")
                return

            # Read validation file
            validations_data: List[Dict[str, Any]] = []
            with open(validation_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    validations_data.append(dict(row))

            if not validations_data:
                ui.label("No validations found").classes("text-gray-500 text-lg italic")
                return

            # Filter state - all statuses selected by default
            all_validation_statuses = ["present", "absent", "uncertain", "conflicting"]
            filter_validations: Dict[str, List[str]] = {
                "value": list(all_validation_statuses)
            }

            # Create filter menu first (before table)
            create_validation_filter_menu(
                all_statuses=all_validation_statuses,
                filter_state=filter_validations,
                on_change=lambda: refresh_table(),
            )

            # Table container
            table_container = ui.column().classes("w-full")

            @ui.refreshable
            def refresh_table():
                """Refresh the table with current filters."""
                table_container.clear()

                # Apply filters
                filtered_data = validations_data.copy()
                if filter_validations["value"]:
                    filtered_data = [
                        row
                        for row in filtered_data
                        if row.get("Validation") in filter_validations["value"]
                    ]

                with table_container:
                    # Show count
                    if filter_validations["value"] != all_validation_statuses:
                        ui.label(
                            f"Showing {len(filtered_data)} of {len(validations_data)} validations"
                        ).classes("text-sm text-gray-600 mb-2")
                    else:
                        ui.label(f"{len(filtered_data)} validations").classes(
                            "text-sm text-gray-600 mb-2"
                        )

                    # Prepare columns for table
                    columns: List[Dict[str, Any]] = [
                        {"name": "actions", "label": "", "field": "actions"},
                        {
                            "name": "FID",
                            "label": "Family ID",
                            "field": "FID",
                            "sortable": True,
                            "align": "left",
                        },
                        {
                            "name": "Variant",
                            "label": "Variant",
                            "field": "Variant",
                            "sortable": True,
                            "align": "left",
                        },
                        {
                            "name": "Sample",
                            "label": "Sample",
                            "field": "Sample",
                            "sortable": True,
                            "align": "left",
                        },
                        {
                            "name": "User",
                            "label": "User",
                            "field": "User",
                            "sortable": True,
                            "align": "left",
                        },
                        {
                            "name": "Inheritance",
                            "label": "Inheritance",
                            "field": "Inheritance",
                            "sortable": True,
                            "align": "left",
                        },
                        {
                            "name": "Validation",
                            "label": "Validation",
                            "field": "Validation",
                            "sortable": True,
                            "align": "left",
                        },
                        {
                            "name": "Timestamp",
                            "label": "Timestamp",
                            "field": "Timestamp",
                            "sortable": True,
                            "align": "left",
                        },
                    ]

                    # Create table
                    validation_table = (
                        ui.table(
                            columns=columns,
                            rows=filtered_data,
                            row_key="Timestamp",
                            pagination={"rowsPerPage": 50},
                        )
                        .classes("w-full")
                        .props("dense flat")
                    )

                    # Add custom slot for view button and validation icons
                    validation_table.add_slot("body", VALIDATION_TABLE_SLOT)

                    # Handle view button click
                    def on_view_variant(e):
                        row_data = e.args
                        family_id = row_data.get("FID", "")
                        variant_str = row_data.get("Variant", "")
                        sample_id = row_data.get("Sample", "")

                        try:
                            parts = variant_str.split(":")
                            if len(parts) == 4:
                                chrom, pos, ref, alt = parts

                                # Find the cohort from family_id
                                cohort_name = None
                                for c_name, cohort in store.cohorts.items():
                                    if family_id in cohort.families:
                                        cohort_name = c_name
                                        break

                                if not cohort_name:
                                    ui.notify(
                                        f"Could not find cohort for family {family_id}",
                                        type="warning",
                                    )
                                    return

                                # Create variant data dict
                                variant_data = dict(row_data)

                                # Callback to update the Validation column in the table
                                def on_save(validation_status: str):
                                    # Reload validation data from file
                                    validations_data.clear()
                                    with open(validation_file, "r") as f:
                                        reader = csv.DictReader(f, delimiter="\t")
                                        for row in reader:
                                            validations_data.append(dict(row))
                                    # Refresh the table display
                                    refresh_table()

                                # Show dialog
                                show_variant_dialog(
                                    cohort_name=cohort_name,
                                    family_id=family_id,
                                    chrom=chrom,
                                    pos=pos,
                                    ref=ref,
                                    alt=alt,
                                    sample=sample_id,
                                    variant_data=variant_data,
                                    on_save_callback=on_save,
                                )
                            else:
                                ui.notify(
                                    "Invalid variant format. Expected chr:pos:ref:alt",
                                    type="warning",
                                )
                        except Exception as ex:
                            ui.notify(f"Error parsing variant: {ex}", type="warning")

                    validation_table.on("view_variant", on_view_variant)

            # Initial render
            refresh_table()

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
