"""Validation statistics page."""

import csv
from collections import Counter
from datetime import datetime
from typing import Dict, List

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.utils.data import get_data_store


@ui.page("/validation/statistics")
def validation_statistics_page() -> None:
    """Render the validation statistics page."""
    create_header()

    try:
        store = get_data_store()
        snv_validation_file = store.data_dir / "validations" / "snvs.tsv"
        sv_validation_file = store.data_dir / "validations" / "svs.tsv"

        with ui.column().classes("w-full px-6 py-6"):
            ui.label("📊 Validation Statistics").classes(
                "text-3xl font-bold text-blue-900 mb-6"
            )

            # Load validation data from both SNV and SV files
            def load_all_validations():
                """Load validations from both SNV and SV files."""
                all_validations = []
                ignored_count = 0

                # Load SNV validations
                if snv_validation_file.exists():
                    with open(snv_validation_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            row_dict = dict(row)
                            row_dict["Type"] = "SNV"
                            if row.get("Ignore", "0") == "1":
                                ignored_count += 1
                            else:
                                all_validations.append(row_dict)

                # Load SV validations
                if sv_validation_file.exists():
                    with open(sv_validation_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            row_dict = dict(row)
                            row_dict["Type"] = "SV"
                            if row.get("Ignore", "0") == "1":
                                ignored_count += 1
                            else:
                                all_validations.append(row_dict)

                return all_validations, ignored_count

            all_validations_data, total_ignored_count = load_all_validations()

            if not all_validations_data:
                ui.label("No validation data available").classes(
                    "text-gray-500 text-lg italic"
                )
            else:
                # Extract all unique users
                all_unique_users = sorted(
                    set(
                        row.get("User", "")
                        for row in all_validations_data
                        if row.get("User")
                    )
                )

                # Filter states
                filter_types: Dict[str, List[str]] = {"value": ["SNV", "SV"]}
                filter_users: Dict[str, List[str]] = {"value": list(all_unique_users)}

                # Create filter UI
                with ui.row().classes("gap-4 mb-6 flex-wrap items-start"):
                    # Type filter
                    with ui.column().classes("gap-1"):
                        ui.label("Filter by Type:").classes(
                            "text-xs font-semibold text-gray-600"
                        )
                        with ui.row().classes("gap-2"):
                            type_snv = ui.checkbox("SNV", value=True).classes("text-sm")
                            type_sv = ui.checkbox("SV", value=True).classes("text-sm")

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

                                    ui.separator()
                                    with ui.row().classes("gap-2"):

                                        def select_all_users():
                                            for cb in user_checkboxes.values():
                                                cb.value = True
                                            filter_users["value"] = list(
                                                all_unique_users
                                            )
                                            refresh_statistics()

                                        def select_no_users():
                                            for cb in user_checkboxes.values():
                                                cb.value = False
                                            filter_users["value"] = []
                                            refresh_statistics()

                                        ui.button(
                                            "All", on_click=select_all_users
                                        ).props("flat dense size=sm")
                                        ui.button(
                                            "None", on_click=select_no_users
                                        ).props("flat dense size=sm")

                # Container for statistics
                stats_container = ui.column().classes("w-full")

                @ui.refreshable
                def refresh_statistics():
                    """Refresh statistics based on filters."""
                    stats_container.clear()

                    # Apply filters
                    filtered_data = all_validations_data.copy()

                    # Filter by type
                    if filter_types["value"]:
                        filtered_data = [
                            row
                            for row in filtered_data
                            if row.get("Type") in filter_types["value"]
                        ]

                    # Filter by user
                    if filter_users["value"]:
                        filtered_data = [
                            row
                            for row in filtered_data
                            if row.get("User") in filter_users["value"]
                        ]

                    validations_data = filtered_data
                    ignored_count = (
                        total_ignored_count  # Note: ignored count is not filtered
                    )

                    with stats_container:
                        # Overall statistics with Type pie chart
                        with ui.card().classes("w-full mb-6 p-4"):
                            ui.label("Overall Statistics").classes(
                                "text-xl font-semibold mb-4 text-blue-800"
                            )
                            with ui.row().classes("gap-8 flex-wrap items-start"):
                                # Statistics numbers
                                with ui.column().classes("gap-4"):
                                    with ui.row().classes("gap-8 flex-wrap"):
                                        with ui.column().classes("gap-1"):
                                            ui.label("Total Validations").classes(
                                                "text-sm text-gray-600"
                                            )
                                            ui.label(
                                                str(len(validations_data))
                                            ).classes(
                                                "text-3xl font-bold text-blue-700"
                                            )

                                        unique_variants = len(
                                            set(
                                                row.get("Variant", "")
                                                for row in validations_data
                                            )
                                        )
                                        with ui.column().classes("gap-1"):
                                            ui.label("Unique Variants").classes(
                                                "text-sm text-gray-600"
                                            )
                                            ui.label(str(unique_variants)).classes(
                                                "text-3xl font-bold text-green-700"
                                            )

                                        unique_families = len(
                                            set(
                                                row.get("FID", "")
                                                for row in validations_data
                                            )
                                        )
                                        with ui.column().classes("gap-1"):
                                            ui.label("Families").classes(
                                                "text-sm text-gray-600"
                                            )
                                            ui.label(str(unique_families)).classes(
                                                "text-3xl font-bold text-purple-700"
                                            )

                                        unique_samples = len(
                                            set(
                                                row.get("Sample", "")
                                                for row in validations_data
                                            )
                                        )
                                        with ui.column().classes("gap-1"):
                                            ui.label("Samples").classes(
                                                "text-sm text-gray-600"
                                            )
                                            ui.label(str(unique_samples)).classes(
                                                "text-3xl font-bold text-orange-700"
                                            )

                                        # Show ignored count if any
                                        if ignored_count > 0:
                                            with ui.column().classes("gap-1"):
                                                ui.label("Ignored").classes(
                                                    "text-sm text-gray-600"
                                                )
                                                ui.label(str(ignored_count)).classes(
                                                    "text-3xl font-bold text-gray-400"
                                                )

                        # Charts in a grid
                        with ui.row().classes("w-full gap-4 flex-wrap"):
                            # Validation Status Chart
                            with ui.card().classes("flex-1 min-w-[400px] p-4"):
                                ui.label("Validation Status").classes(
                                    "text-lg font-semibold mb-2 text-gray-800"
                                )
                                validation_counts = Counter(
                                    row.get("Validation", "Unknown")
                                    for row in validations_data
                                )
                                ui.echart(
                                    {
                                        "tooltip": {"trigger": "item"},
                                        "series": [
                                            {
                                                "type": "pie",
                                                "radius": "70%",
                                                "data": [
                                                    {
                                                        "name": status,
                                                        "value": count,
                                                        "itemStyle": {
                                                            "color": {
                                                                "present": "#22c55e",
                                                                "in phase MNV": "#16a34a",
                                                                "absent": "#ef4444",
                                                                "uncertain": "#f59e0b",
                                                                "different": "#fb923c",
                                                                "conflicting": "#fbbf24",
                                                            }.get(status, "#6b7280")
                                                        },
                                                    }
                                                    for status, count in validation_counts.items()
                                                ],
                                                "label": {"formatter": "{b}: {c}"},
                                            }
                                        ],
                                    }
                                ).classes("w-full h-64")

                            # Type Distribution Chart
                            with ui.card().classes("flex-1 min-w-[400px] p-4"):
                                ui.label("Type Distribution").classes(
                                    "text-lg font-semibold mb-2 text-gray-800"
                                )
                                type_counts = Counter(
                                    row.get("Type", "Unknown")
                                    for row in validations_data
                                )
                                ui.echart(
                                    {
                                        "tooltip": {"trigger": "item"},
                                        "series": [
                                            {
                                                "type": "pie",
                                                "radius": "70%",
                                                "data": [
                                                    {
                                                        "name": vtype,
                                                        "value": count,
                                                        "itemStyle": {
                                                            "color": {
                                                                "SNV": "#3b82f6",
                                                                "SV": "#8b5cf6",
                                                            }.get(vtype, "#6b7280")
                                                        },
                                                    }
                                                    for vtype, count in type_counts.items()
                                                ],
                                                "label": {
                                                    "formatter": "{b}: {c} ({d}%)"
                                                },
                                            }
                                        ],
                                    }
                                ).classes("w-full h-64")

                        # User Activity and Timeline
                        with ui.row().classes("w-full gap-4 flex-wrap mt-4"):
                            # Validations by User
                            with ui.card().classes("flex-1 min-w-[400px] p-4"):
                                ui.label("Validations by User").classes(
                                    "text-lg font-semibold mb-2 text-gray-800"
                                )
                                user_counts = Counter(
                                    row.get("User", "Unknown")
                                    for row in validations_data
                                )
                                ui.echart(
                                    {
                                        "tooltip": {},
                                        "xAxis": {
                                            "type": "category",
                                            "data": list(user_counts.keys()),
                                            "axisLabel": {"rotate": -45},
                                        },
                                        "yAxis": {"type": "value", "name": "Count"},
                                        "series": [
                                            {
                                                "type": "bar",
                                                "data": list(user_counts.values()),
                                                "itemStyle": {"color": "#8b5cf6"},
                                            }
                                        ],
                                    }
                                ).classes("w-full h-64")

                            # Inheritance Patterns
                            with ui.card().classes("flex-1 min-w-[400px] p-4"):
                                ui.label("Inheritance Patterns").classes(
                                    "text-lg font-semibold mb-2 text-gray-800"
                                )
                                inheritance_counts = Counter(
                                    row.get("Inheritance", "Unknown")
                                    for row in validations_data
                                )
                                ui.echart(
                                    {
                                        "tooltip": {},
                                        "xAxis": {
                                            "type": "category",
                                            "data": list(inheritance_counts.keys()),
                                            "axisLabel": {"rotate": -45},
                                        },
                                        "yAxis": {"type": "value", "name": "Count"},
                                        "series": [
                                            {
                                                "type": "bar",
                                                "data": list(
                                                    inheritance_counts.values()
                                                ),
                                                "itemStyle": {"color": "#3b82f6"},
                                            }
                                        ],
                                    }
                                ).classes("w-full h-64")

                        # Timeline chart (if timestamps are available)
                        with ui.card().classes("w-full p-4 mt-4"):
                            ui.label("Validation Timeline").classes(
                                "text-lg font-semibold mb-2 text-gray-800"
                            )
                            try:
                                # Parse timestamps and group by date
                                date_counts: Counter = Counter()
                                for row in validations_data:
                                    timestamp_str = row.get("Timestamp", "")
                                    if timestamp_str:
                                        try:
                                            # Parse ISO format timestamp
                                            dt = datetime.fromisoformat(
                                                timestamp_str.replace("Z", "+00:00")
                                            )
                                            date_key = dt.strftime("%Y-%m-%d")
                                            date_counts[date_key] += 1
                                        except Exception:
                                            pass

                                if date_counts:
                                    sorted_dates = sorted(date_counts.keys())
                                    ui.echart(
                                        {
                                            "tooltip": {"trigger": "axis"},
                                            "xAxis": {
                                                "type": "category",
                                                "data": sorted_dates,
                                                "name": "Date",
                                            },
                                            "yAxis": {"type": "value", "name": "Count"},
                                            "series": [
                                                {
                                                    "type": "line",
                                                    "data": [
                                                        date_counts[date]
                                                        for date in sorted_dates
                                                    ],
                                                    "smooth": True,
                                                    "itemStyle": {"color": "#10b981"},
                                                    "lineStyle": {"color": "#10b981"},
                                                }
                                            ],
                                        }
                                    ).classes("w-full h-64")
                                else:
                                    ui.label("No timeline data available").classes(
                                        "text-gray-500 italic"
                                    )
                            except Exception as e:
                                ui.label(f"Could not parse timeline: {e}").classes(
                                    "text-gray-500 italic"
                                )

                # Setup filter handlers
                def on_type_change():
                    filter_types["value"] = []
                    if type_snv.value:
                        filter_types["value"].append("SNV")
                    if type_sv.value:
                        filter_types["value"].append("SV")
                    refresh_statistics()

                type_snv.on_value_change(lambda: on_type_change())
                type_sv.on_value_change(lambda: on_type_change())

                for user, cb in user_checkboxes.items():

                    def make_user_handler(u):
                        def handler():
                            filter_users["value"] = [
                                user for user, cb in user_checkboxes.items() if cb.value
                            ]
                            refresh_statistics()

                        return handler

                    cb.on_value_change(make_user_handler(user))

                # Initial render
                refresh_statistics()

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
