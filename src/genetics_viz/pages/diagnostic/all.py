"""Diagnostic all page - displays all diagnostics from diagnostics/snvs.tsv and svs.tsv."""

import asyncio
import csv
from datetime import datetime
from typing import Any, Dict, List, Optional

from nicegui import ui

from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.header import create_header
from genetics_viz.components.icons import get_diagnostic_icon
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.auth import check_auth
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.diagnostic_badges import build_diagnostic_badge
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.vep import format_consequence_display, get_consequence_color
from genetics_viz.utils.wisecondorx import WISECONDORX_CONFIG


@ui.page("/diagnostic/all")
async def diagnostic_all_page() -> None:
    """Render all diagnostics from diagnostics/snvs.tsv and svs.tsv."""
    if redirect := check_auth():
        return redirect
    create_header()

    try:
        store = get_data_store()
        snv_diagnostic_file = store.data_dir / "diagnostics" / "snvs.tsv"
        sv_diagnostic_file = store.data_dir / "diagnostics" / "svs.tsv"

        with ui.column().classes("w-full px-6 py-6"):
            # Title
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label("🏥 All Diagnostics").classes(
                    "text-3xl font-bold text-blue-900"
                )

            if not snv_diagnostic_file.exists() and not sv_diagnostic_file.exists():
                ui.label("No diagnostics found").classes("text-gray-500 text-lg italic")
                return

            # Function to load and aggregate diagnostic data
            def load_and_aggregate_diagnostics() -> List[Dict[str, Any]]:
                """Load raw diagnostics from both files and aggregate by (Type, FID, Variant, Sample)."""
                raw_diagnostics: List[Dict[str, Any]] = []

                # Load SNV diagnostics
                if snv_diagnostic_file.exists():
                    with open(snv_diagnostic_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            row_dict = dict(row)
                            row_dict["Type"] = "SNV"
                            raw_diagnostics.append(row_dict)

                # Load SV diagnostics
                if sv_diagnostic_file.exists():
                    with open(sv_diagnostic_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            row_dict = dict(row)
                            row_dict["Type"] = "SV"
                            raw_diagnostics.append(row_dict)

                if not raw_diagnostics:
                    return []

                # Aggregate diagnostics by (Type, FID, Variant, Sample)
                aggregated: Dict[tuple, Dict[str, Any]] = {}

                for row in raw_diagnostics:
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
                            "Gene": row.get("Gene", ""),
                            "Impact": row.get("Impact", ""),
                            "Sample": row.get("Sample", ""),
                            "users": set(),
                            "diagnostics": [],
                            "timestamps": [],
                        }

                    agg = aggregated[key]
                    agg["users"].add(row.get("User", ""))
                    agg["diagnostics"].append(
                        (
                            row.get("Diagnostic", ""),
                            row.get("User", ""),
                            row.get("Timestamp", ""),
                            row.get("Comment", ""),
                            row.get("Ignore", "0"),
                        )
                    )
                    agg["timestamps"].append(row.get("Timestamp", ""))

                # Convert to list with aggregated values
                diagnostics_data: List[Dict[str, Any]] = []
                for _key, data in aggregated.items():
                    # Filter non-ignored entries
                    non_ignored = [d for d in data["diagnostics"] if d[4] != "1"]

                    if not non_ignored:
                        continue

                    diag_values = [d[0] for d in non_ignored]
                    unique_diags = set(diag_values)

                    # Determine final diagnostic status
                    if len(unique_diags) > 1:
                        final_diagnostic = "conflicting"
                    else:
                        final_diagnostic = diag_values[0]

                    # Build badge tuples: (diagnostic_value, user, timestamp, comment)
                    badge_data = [(d[0], d[1], d[2], d[3]) for d in non_ignored]

                    diagnostics_data.append(
                        {
                            "Type": data["Type"],
                            "FID": data["FID"],
                            "Variant": data["Variant"],
                            "Gene": data["Gene"],
                            "Impact": data["Impact"],
                            "Sample": data["Sample"],
                            "User": ", ".join(sorted(data["users"])),
                            "Diagnostic": final_diagnostic,
                            "Diagnostic_badge": build_diagnostic_badge(
                                final_diagnostic, badge_data
                            ),
                            "Timestamp": max(data["timestamps"]),
                        }
                    )

                # Build gene badges and impact badges for each row
                gene_scorer = get_gene_scorer()
                for row in diagnostics_data:
                    # Gene badges
                    gene_str = row.get("Gene", "")
                    if gene_str:
                        gene_badges = []
                        # SNVs: single symbol; SVs: comma-separated
                        for part in str(gene_str).split(","):
                            symbol = part.strip()
                            if symbol:
                                score, _ = gene_scorer.get_gene_score_and_sets(symbol)
                                color = gene_scorer.get_gene_color(symbol)
                                tooltip = gene_scorer.get_gene_tooltip(symbol)
                                gene_badges.append(
                                    {
                                        "label": symbol,
                                        "color": color,
                                        "tooltip": tooltip,
                                        "score": score,
                                    }
                                )
                        # Sort by score descending
                        gene_badges.sort(key=lambda x: x["score"], reverse=True)
                        # Truncate to 6 + "+X genes"
                        total = len(gene_badges)
                        if total > 6:
                            gene_badges = gene_badges[:6]
                            gene_badges.append(
                                {
                                    "label": f"+{total - 6} genes",
                                    "color": "#9e9e9e",
                                    "tooltip": f"{total - 6} more genes",
                                }
                            )
                        row["GeneBadges"] = gene_badges
                    else:
                        row["GeneBadges"] = []

                    # Impact badges
                    impact_str = row.get("Impact", "")
                    row_type = row.get("Type", "")
                    if impact_str:
                        if row_type == "SV":
                            # SV: GAIN/LOSS with colors from WisecondorX config
                            upper = impact_str.upper()
                            if "GAIN" in upper:
                                row["ImpactBadges"] = [
                                    {
                                        "label": "GAIN",
                                        "color": WISECONDORX_CONFIG["robust_gain"][
                                            "color"
                                        ],
                                    }
                                ]
                            elif "LOSS" in upper:
                                row["ImpactBadges"] = [
                                    {
                                        "label": "LOSS",
                                        "color": WISECONDORX_CONFIG["robust_loss"][
                                            "color"
                                        ],
                                    }
                                ]
                            else:
                                row["ImpactBadges"] = [
                                    {"label": impact_str, "color": "#6b7280"}
                                ]
                        else:
                            # SNV: VEP consequence badges
                            # Split by both ',' and '&' to handle aggregated values
                            impact_badges: List[Dict[str, str]] = []
                            seen: set = set()
                            for part in str(impact_str).split(","):
                                for cons in part.split("&"):
                                    cons = cons.strip()
                                    if cons:
                                        label = format_consequence_display(cons)
                                        color = get_consequence_color(cons)
                                        key = (label, color)
                                        if key not in seen:
                                            seen.add(key)
                                            impact_badges.append(
                                                {"label": label, "color": color}
                                            )
                            row["ImpactBadges"] = impact_badges
                    else:
                        row["ImpactBadges"] = []

                return diagnostics_data

            # Load initial data (offloaded to thread)
            diagnostics_data = await asyncio.to_thread(load_and_aggregate_diagnostics)

            if not diagnostics_data:
                ui.label("No diagnostics found").classes("text-gray-500 text-lg italic")
                return

            # Extract all unique users
            all_users = sorted(
                set(row["User"] for row in diagnostics_data if row["User"])
            )
            all_unique_users: List[str] = sorted(
                {
                    user.strip()
                    for user_str in all_users
                    for user in user_str.split(", ")
                }
            )

            # Filter state
            all_diagnostic_statuses = [
                "pathogenic",
                "uncertain",
                "benign",
                "conflicting",
            ]
            filter_diagnostics: Dict[str, List[str]] = {
                "value": list(all_diagnostic_statuses)
            }
            filter_types: Dict[str, List[str]] = {"value": ["SNV", "SV"]}
            filter_users: Dict[str, List[str]] = {"value": list(all_unique_users)}

            # Date filter state
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

                # Diagnostic status filter
                create_validation_filter_menu(
                    all_statuses=all_diagnostic_statuses,
                    filter_state=filter_diagnostics,
                    on_change=lambda: refresh_table(),
                    label="Filter by Diagnostic",
                    icon_fn=get_diagnostic_icon,
                )

                # User filter
                with ui.column().classes("gap-1"):
                    ui.label("Filter by User:").classes(
                        "text-xs font-semibold text-gray-600"
                    )
                    with ui.button(icon="person", color="blue").props("flat dense"):
                        with ui.menu():
                            with ui.column().classes("p-2 gap-1"):
                                user_checkboxes: Dict[str, Any] = {}
                                for user in all_unique_users:
                                    cb = ui.checkbox(user, value=True).classes(
                                        "text-sm"
                                    )
                                    user_checkboxes[user] = cb

                                    def make_user_handler(u):
                                        def handler():
                                            filter_users["value"] = [
                                                usr
                                                for usr, ucb in user_checkboxes.items()
                                                if ucb.value
                                            ]
                                            refresh_table()

                                        return handler

                                    cb.on_value_change(make_user_handler(user))

                                ui.separator()
                                with ui.row().classes("gap-2"):

                                    def select_all_users():
                                        for ucb in user_checkboxes.values():
                                            ucb.value = True
                                        filter_users["value"] = list(all_unique_users)
                                        refresh_table()

                                    def select_no_users():
                                        for ucb in user_checkboxes.values():
                                            ucb.value = False
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

            @ui.refreshable
            def refresh_table():
                """Refresh the table with current filters."""
                table_container.clear()

                # Apply filters
                filtered_data = diagnostics_data.copy()

                # Filter by type
                if filter_types["value"]:
                    filtered_data = [
                        row
                        for row in filtered_data
                        if row.get("Type") in filter_types["value"]
                    ]

                # Filter by diagnostic status
                if filter_diagnostics["value"]:
                    filtered_data = [
                        row
                        for row in filtered_data
                        if row.get("Diagnostic") in filter_diagnostics["value"]
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

                    def parse_timestamp(ts_str: str):
                        try:
                            dt = datetime.fromisoformat(ts_str.replace("Z", "+00:00"))
                            return dt.date()
                        except Exception:
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
                    if filter_diagnostics["value"] != all_diagnostic_statuses:
                        active_filters.append(
                            f"Status: {', '.join(filter_diagnostics['value'])}"
                        )
                    if filter_users["value"] != list(all_unique_users):
                        active_filters.append(
                            f"Users: {len(filter_users['value'])} selected"
                        )
                    if filter_date_mode["value"] != "all":
                        active_filters.append(f"Date: {date_mode_select.value}")

                    if active_filters:
                        ui.label(
                            f"Showing {len(filtered_data)} of {len(diagnostics_data)}"
                            f" diagnostics ({'; '.join(active_filters)})"
                        ).classes("text-sm text-gray-600 mb-2")
                    else:
                        ui.label(f"{len(filtered_data)} diagnostics").classes(
                            "text-sm text-gray-600 mb-2"
                        )

                    # Prepare columns for table
                    columns: List[Dict[str, Any]] = [
                        {"id": "Type", "header": "Type", "sortable": True},
                        {"id": "FID", "header": "Family ID", "sortable": True},
                        {
                            "id": "Variant",
                            "header": "Variant",
                            "sortable": True,
                        },
                        {
                            "id": "Gene",
                            "header": "Gene",
                            "sortable": True,
                            "cellType": "gene_badge",
                            "badgesField": "GeneBadges",
                        },
                        {
                            "id": "Impact",
                            "header": "Impact",
                            "sortable": True,
                            "cellType": "badge_list",
                            "badgesField": "ImpactBadges",
                        },
                        {"id": "Sample", "header": "Sample", "sortable": True},
                        {"id": "User", "header": "User", "sortable": True},
                        {
                            "id": "Diagnostic",
                            "header": "Diagnostic",
                            "cellType": "diagnostic",
                            "sortable": True,
                        },
                        {
                            "id": "Timestamp",
                            "header": "Timestamp",
                            "sortable": True,
                        },
                    ]

                    DataTable(
                        columns=columns,
                        rows=filtered_data,
                        row_key="Timestamp",
                        pagination={"rowsPerPage": 50},
                    )

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
