"""Validation all page - displays all validations from validations/snvs.tsv and svs.tsv."""

import csv
from datetime import datetime
from typing import Any, Dict, List, Optional

from nicegui import app as nicegui_app
from nicegui import ui

from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.header import create_header
from genetics_viz.components.sv_dialog import show_sv_dialog
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
        snv_validation_file = store.data_dir / "validations" / "snvs.tsv"
        sv_validation_file = store.data_dir / "validations" / "svs.tsv"

        # Serve data files for IGV.js
        nicegui_app.add_static_files("/data", str(store.data_dir))

        with ui.column().classes("w-full px-6 py-6"):
            # Title
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label("📋 All Validations").classes(
                    "text-3xl font-bold text-blue-900"
                )

            if not snv_validation_file.exists() and not sv_validation_file.exists():
                ui.label("No validations found").classes("text-gray-500 text-lg italic")
                return

            # Function to load and aggregate validation data
            def load_and_aggregate_validations() -> List[Dict[str, Any]]:
                """Load raw validations from both SNV and SV files and aggregate by (Type, FID, Variant, Sample)."""
                raw_validations: List[Dict[str, Any]] = []

                # Load SNV validations
                if snv_validation_file.exists():
                    with open(snv_validation_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            row_dict = dict(row)
                            row_dict["Type"] = "SNV"
                            raw_validations.append(row_dict)

                # Load SV validations
                if sv_validation_file.exists():
                    with open(sv_validation_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            row_dict = dict(row)
                            row_dict["Type"] = "SV"
                            # For SVs, format variant with curated coordinates if available
                            curated_start = row.get("CuratedStart", "")
                            curated_end = row.get("CuratedEnd", "")

                            # Parse original variant (chr:start-end:type)
                            variant_str = row.get("Variant", "")
                            if curated_start and curated_end and ":" in variant_str:
                                parts = variant_str.split(":")
                                if len(parts) >= 3:
                                    chrom = parts[0]
                                    sv_type = parts[2] if len(parts) == 3 else parts[3]
                                    # Store original for lookup
                                    row_dict["OriginalVariant"] = variant_str
                                    # Display curated in table
                                    row_dict["Variant"] = (
                                        f"{chrom}:{curated_start}-{curated_end}:{sv_type}"
                                    )
                                    row_dict["IsCurated"] = True
                                else:
                                    row_dict["OriginalVariant"] = variant_str
                                    row_dict["IsCurated"] = False
                            else:
                                row_dict["OriginalVariant"] = variant_str
                                row_dict["IsCurated"] = False

                            raw_validations.append(row_dict)

                if not raw_validations:
                    return []

                # Aggregate validations by (Type, FID, Variant, Sample)
                aggregated: Dict[tuple, Dict[str, Any]] = {}

                for row in raw_validations:
                    key = (
                        row.get("Type", ""),
                        row.get("FID", ""),
                        row.get("Variant", ""),
                        row.get("Sample", ""),
                    )

                    if key not in aggregated:
                        aggregated[key] = {
                            "Type": row.get("Type", ""),
                            "FID": row.get("FID", ""),
                            "Variant": row.get("Variant", ""),
                            "Sample": row.get("Sample", ""),
                            "OriginalVariant": row.get(
                                "OriginalVariant", row.get("Variant", "")
                            ),
                            "IsCurated": row.get("IsCurated", False),
                            "users": set(),
                            "validations": [],
                            "inheritances": [],
                            "timestamps": [],
                            "ignored": [],
                        }

                    # Collect data
                    aggregated[key]["users"].add(row.get("User", ""))
                    aggregated[key]["validations"].append(row.get("Validation", ""))
                    aggregated[key]["inheritances"].append(row.get("Inheritance", ""))
                    aggregated[key]["timestamps"].append(row.get("Timestamp", ""))
                    aggregated[key]["ignored"].append(row.get("Ignore", "0"))

                # Convert to list with aggregated values
                validations_data: List[Dict[str, Any]] = []
                for key, data in aggregated.items():
                    # Filter non-ignored validations for status computation
                    non_ignored_vals = [
                        (v, i)
                        for v, i, ig in zip(
                            data["validations"], data["inheritances"], data["ignored"]
                        )
                        if ig != "1"
                    ]

                    if not non_ignored_vals:
                        # All are ignored, skip
                        continue

                    validation_statuses = [v[0] for v in non_ignored_vals]
                    # Normalize "in phase MNV" to "present" for conflict detection
                    normalized_statuses = [
                        "present" if s == "in phase MNV" else s
                        for s in validation_statuses
                    ]
                    unique_validations = set(normalized_statuses)

                    # Determine final validation status
                    if len(unique_validations) > 1:
                        final_validation = "conflicting"
                        final_inheritance = ""
                    elif "present" in unique_validations:
                        # Check if any is specifically "in phase MNV"
                        if "in phase MNV" in validation_statuses:
                            final_validation = "in phase MNV"
                        else:
                            final_validation = "present"
                        # Check inheritance - prioritize de novo, then homozygous
                        inheritances = [
                            v[1]
                            for v in non_ignored_vals
                            if v[0] in ("present", "in phase MNV")
                        ]
                        if "de novo" in inheritances:
                            final_inheritance = "de novo"
                        elif "homozygous" in inheritances:
                            final_inheritance = "homozygous"
                        else:
                            final_inheritance = ""
                    elif "absent" in unique_validations:
                        final_validation = "absent"
                        final_inheritance = ""
                    else:
                        final_validation = "uncertain"
                        final_inheritance = ""

                    validations_data.append(
                        {
                            "Type": data["Type"],
                            "FID": data["FID"],
                            "Variant": data["Variant"],
                            "OriginalVariant": data["OriginalVariant"],
                            "IsCurated": data["IsCurated"],
                            "Sample": data["Sample"],
                            "User": ", ".join(sorted(data["users"])),
                            "Inheritance": final_inheritance,
                            "Validation": final_validation,
                            "Timestamp": max(
                                data["timestamps"]
                            ),  # Most recent timestamp
                        }
                    )

                return validations_data

            # Load initial data
            validations_data = load_and_aggregate_validations()

            if not validations_data:
                ui.label("No validations found").classes("text-gray-500 text-lg italic")
                return

            # Extract all unique users
            all_users = sorted(
                set(row["User"] for row in validations_data if row["User"])
            )
            all_unique_users = set()
            for user_str in all_users:
                for user in user_str.split(", "):
                    all_unique_users.add(user.strip())
            all_unique_users = sorted(all_unique_users)

            # Filter state - all statuses selected by default
            all_validation_statuses = ["present", "absent", "uncertain", "conflicting"]
            filter_validations: Dict[str, List[str]] = {
                "value": list(all_validation_statuses)
            }

            # Filter by type - default both
            filter_types: Dict[str, List[str]] = {"value": ["SNV", "SV"]}

            # Filter by user - default all
            filter_users: Dict[str, List[str]] = {"value": list(all_unique_users)}

            # Filter by date - default show all
            filter_date_mode: Dict[str, str] = {"value": "all"}
            filter_date_before: Dict[str, Optional[str]] = {"value": None}
            filter_date_after: Dict[str, Optional[str]] = {"value": None}
            filter_date_between_start: Dict[str, Optional[str]] = {"value": None}
            filter_date_between_end: Dict[str, Optional[str]] = {"value": None}

            # Create filter UI
            with ui.row().classes("gap-4 mb-4 flex-wrap items-start"):
                # Type filter
                with ui.column().classes("gap-1"):
                    ui.label("Filter by Type:").classes(
                        "text-xs font-semibold text-gray-600"
                    )
                    with ui.row().classes("gap-2"):
                        type_snv = ui.checkbox("SNV", value=True).classes("text-sm")
                        type_sv = ui.checkbox("SV", value=True).classes("text-sm")

                        def on_type_change():
                            filter_types["value"] = []
                            if type_snv.value:
                                filter_types["value"].append("SNV")
                            if type_sv.value:
                                filter_types["value"].append("SV")
                            refresh_table()

                        type_snv.on_value_change(lambda: on_type_change())
                        type_sv.on_value_change(lambda: on_type_change())

                # Validation status filter
                create_validation_filter_menu(
                    all_statuses=all_validation_statuses,
                    filter_state=filter_validations,
                    on_change=lambda: refresh_table(),
                )

                # User filter
                with ui.column().classes("gap-1"):
                    ui.label("Filter by User:").classes(
                        "text-xs font-semibold text-gray-600"
                    )
                    with ui.button(icon="person", color="blue").props("flat dense"):
                        with ui.menu() as user_menu:
                            with ui.column().classes("p-2 gap-1"):
                                user_checkboxes = {}
                                for user in all_unique_users:
                                    cb = ui.checkbox(user, value=True).classes(
                                        "text-sm"
                                    )
                                    user_checkboxes[user] = cb

                                    def make_user_handler(u):
                                        def handler():
                                            filter_users["value"] = [
                                                user
                                                for user, cb in user_checkboxes.items()
                                                if cb.value
                                            ]
                                            refresh_table()

                                        return handler

                                    cb.on_value_change(make_user_handler(user))

                                ui.separator()
                                with ui.row().classes("gap-2"):

                                    def select_all_users():
                                        for cb in user_checkboxes.values():
                                            cb.value = True
                                        filter_users["value"] = list(all_unique_users)
                                        refresh_table()

                                    def select_no_users():
                                        for cb in user_checkboxes.values():
                                            cb.value = False
                                        filter_users["value"] = []
                                        refresh_table()

                                    ui.button("All", on_click=select_all_users).props(
                                        "flat dense size=sm"
                                    )
                                    ui.button("None", on_click=select_no_users).props(
                                        "flat dense size=sm"
                                    )

                # Date filter
                with ui.column().classes("gap-1"):
                    ui.label("Filter by Date:").classes(
                        "text-xs font-semibold text-gray-600"
                    )
                    date_mode_select = (
                        ui.select(
                            options=[
                                "Show All",
                                "Before Date",
                                "After Date",
                                "Between Dates",
                            ],
                            value="Show All",
                        )
                        .classes("w-48")
                        .props("dense")
                    )

                    # Date input containers
                    date_before_container = ui.column().classes("gap-1")
                    date_after_container = ui.column().classes("gap-1")
                    date_between_container = ui.column().classes("gap-1")

                    with date_before_container:
                        date_before_input = ui.date().classes("w-48")

                    with date_after_container:
                        date_after_input = ui.date().classes("w-48")

                    with date_between_container:
                        ui.label("Start:").classes("text-xs")
                        date_between_start_input = ui.date().classes("w-48")
                        ui.label("End:").classes("text-xs")
                        date_between_end_input = ui.date().classes("w-48")

                    # Hide all date inputs initially
                    date_before_container.set_visibility(False)
                    date_after_container.set_visibility(False)
                    date_between_container.set_visibility(False)

                    def on_date_mode_change():
                        mode = date_mode_select.value
                        date_before_container.set_visibility(mode == "Before Date")
                        date_after_container.set_visibility(mode == "After Date")
                        date_between_container.set_visibility(mode == "Between Dates")

                        if mode == "Show All":
                            filter_date_mode["value"] = "all"
                        elif mode == "Before Date":
                            filter_date_mode["value"] = "before"
                        elif mode == "After Date":
                            filter_date_mode["value"] = "after"
                        elif mode == "Between Dates":
                            filter_date_mode["value"] = "between"

                        refresh_table()

                    date_mode_select.on_value_change(lambda: on_date_mode_change())

                    def on_date_input_change():
                        filter_date_before["value"] = date_before_input.value
                        filter_date_after["value"] = date_after_input.value
                        filter_date_between_start["value"] = (
                            date_between_start_input.value
                        )
                        filter_date_between_end["value"] = date_between_end_input.value
                        refresh_table()

                    date_before_input.on_value_change(lambda: on_date_input_change())
                    date_after_input.on_value_change(lambda: on_date_input_change())
                    date_between_start_input.on_value_change(
                        lambda: on_date_input_change()
                    )
                    date_between_end_input.on_value_change(
                        lambda: on_date_input_change()
                    )

            # Table container
            table_container = ui.column().classes("w-full")

            # Capture the client context for use in callbacks
            from nicegui import context

            page_client = context.client

            @ui.refreshable
            def refresh_table():
                """Refresh the table with current filters."""
                table_container.clear()

                # Apply filters
                filtered_data = validations_data.copy()

                # Filter by type
                if filter_types["value"]:
                    filtered_data = [
                        row
                        for row in filtered_data
                        if row.get("Type") in filter_types["value"]
                    ]

                # Filter by validation status
                if filter_validations["value"]:
                    filtered_data = [
                        row
                        for row in filtered_data
                        if row.get("Validation") in filter_validations["value"]
                    ]

                # Filter by user
                if filter_users["value"]:
                    filtered_data = [
                        row
                        for row in filtered_data
                        if any(
                            user in row.get("User", "")
                            for user in filter_users["value"]
                        )
                    ]

                # Filter by date
                if filter_date_mode["value"] != "all":

                    def parse_timestamp(ts_str: str) -> Optional[datetime]:
                        """Parse timestamp string to datetime object (ignoring time)."""
                        try:
                            # Parse timestamp and extract date only
                            dt = datetime.fromisoformat(ts_str.replace("Z", "+00:00"))
                            return dt.date()
                        except:
                            return None

                    if (
                        filter_date_mode["value"] == "before"
                        and filter_date_before["value"]
                    ):
                        before_date = datetime.strptime(
                            filter_date_before["value"], "%Y-%m-%d"
                        ).date()
                        filtered_data = [
                            row
                            for row in filtered_data
                            if parse_timestamp(row.get("Timestamp", ""))
                            and parse_timestamp(row.get("Timestamp", "")) < before_date
                        ]
                    elif (
                        filter_date_mode["value"] == "after"
                        and filter_date_after["value"]
                    ):
                        after_date = datetime.strptime(
                            filter_date_after["value"], "%Y-%m-%d"
                        ).date()
                        filtered_data = [
                            row
                            for row in filtered_data
                            if parse_timestamp(row.get("Timestamp", ""))
                            and parse_timestamp(row.get("Timestamp", "")) > after_date
                        ]
                    elif (
                        filter_date_mode["value"] == "between"
                        and filter_date_between_start["value"]
                        and filter_date_between_end["value"]
                    ):
                        start_date = datetime.strptime(
                            filter_date_between_start["value"], "%Y-%m-%d"
                        ).date()
                        end_date = datetime.strptime(
                            filter_date_between_end["value"], "%Y-%m-%d"
                        ).date()
                        filtered_data = [
                            row
                            for row in filtered_data
                            if parse_timestamp(row.get("Timestamp", ""))
                            and start_date
                            <= parse_timestamp(row.get("Timestamp", ""))
                            <= end_date
                        ]

                with table_container:
                    # Show count
                    active_filters = []
                    if filter_types["value"] != ["SNV", "SV"]:
                        active_filters.append(
                            f"Type: {', '.join(filter_types['value'])}"
                        )
                    if filter_validations["value"] != all_validation_statuses:
                        active_filters.append(
                            f"Status: {', '.join(filter_validations['value'])}"
                        )
                    if filter_users["value"] != list(all_unique_users):
                        active_filters.append(
                            f"Users: {len(filter_users['value'])} selected"
                        )
                    if filter_date_mode["value"] != "all":
                        active_filters.append(f"Date: {date_mode_select.value}")

                    if active_filters:
                        ui.label(
                            f"Showing {len(filtered_data)} of {len(validations_data)} validations ({'; '.join(active_filters)})"
                        ).classes("text-sm text-gray-600 mb-2")
                    else:
                        ui.label(f"{len(filtered_data)} validations").classes(
                            "text-sm text-gray-600 mb-2"
                        )

                    # Prepare columns for table
                    columns: List[Dict[str, Any]] = [
                        {"name": "actions", "label": "", "field": "actions"},
                        {
                            "name": "Type",
                            "label": "Type",
                            "field": "Type",
                            "sortable": True,
                            "align": "left",
                        },
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
                        variant_type = row_data.get("Type", "SNV")
                        family_id = row_data.get("FID", "")
                        variant_str = row_data.get("Variant", "")
                        original_variant_str = row_data.get(
                            "OriginalVariant", variant_str
                        )
                        sample_id = row_data.get("Sample", "")
                        is_curated = row_data.get("IsCurated", False)

                        try:
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

                            # Callback to update the Validation column in the table
                            def on_save(validation_status: str):
                                # Reload and re-aggregate validation data from file
                                nonlocal validations_data
                                validations_data = load_and_aggregate_validations()
                                # Refresh the table display using the captured client context
                                with page_client:
                                    ui.timer(0.1, refresh_table, once=True)

                            if variant_type == "SNV":
                                parts = variant_str.split(":")
                                if len(parts) == 4:
                                    chrom, pos, ref, alt = parts

                                    # Create variant data dict
                                    variant_data = dict(row_data)

                                    # Show SNV dialog
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
                                        "Invalid SNV format. Expected chr:pos:ref:alt",
                                        type="warning",
                                    )
                            elif variant_type == "SV":
                                # Parse original SV variant (chr:start-end:type)
                                parts = original_variant_str.split(":")
                                if len(parts) >= 3:
                                    chrom = parts[0]
                                    start_end = parts[1].split("-")
                                    if len(start_end) == 2:
                                        start, end = start_end

                                        # Create SV data dict
                                        sv_data = dict(row_data)

                                        # Add curated coordinates if available
                                        if is_curated:
                                            curated_parts = variant_str.split(":")
                                            if len(curated_parts) >= 2:
                                                curated_start_end = curated_parts[
                                                    1
                                                ].split("-")
                                                if len(curated_start_end) == 2:
                                                    sv_data["CuratedStart"] = (
                                                        curated_start_end[0]
                                                    )
                                                    sv_data["CuratedEnd"] = (
                                                        curated_start_end[1]
                                                    )

                                        # Show SV dialog with original coordinates
                                        show_sv_dialog(
                                            cohort_name=cohort_name,
                                            family_id=family_id,
                                            chrom=chrom,
                                            start=start,
                                            end=end,
                                            sample=sample_id,
                                            sv_data=sv_data,
                                            on_validation_saved=on_save,
                                        )
                                    else:
                                        ui.notify(
                                            "Invalid SV format. Expected chr:start-end:type",
                                            type="warning",
                                        )
                                else:
                                    ui.notify(
                                        "Invalid SV format. Expected chr:start-end:type",
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
