"""Validation statistics page."""

import csv
from collections import Counter
from datetime import datetime

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.utils.data import get_data_store


@ui.page("/validation/statistics")
def validation_statistics_page() -> None:
    """Render the validation statistics page."""
    create_header()

    try:
        store = get_data_store()
        validation_file = store.data_dir / "validations" / "snvs.tsv"

        with ui.column().classes("w-full px-6 py-6"):
            ui.label("📊 Validation Statistics").classes(
                "text-3xl font-bold text-blue-900 mb-6"
            )

            # Load validation data
            validations_data = []
            ignored_count = 0
            if validation_file.exists():
                with open(validation_file, "r") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    for row in reader:
                        # Separate ignored validations
                        if row.get("Ignore", "0") == "1":
                            ignored_count += 1
                        else:
                            validations_data.append(row)

            if not validations_data:
                ui.label("No validation data available").classes(
                    "text-gray-500 text-lg italic"
                )
            else:
                # Overall statistics
                with ui.card().classes("w-full mb-6 p-4"):
                    ui.label("Overall Statistics").classes(
                        "text-xl font-semibold mb-4 text-blue-800"
                    )
                    with ui.row().classes("gap-8 flex-wrap"):
                        with ui.column().classes("gap-1"):
                            ui.label("Total Validations").classes(
                                "text-sm text-gray-600"
                            )
                            ui.label(str(len(validations_data))).classes(
                                "text-3xl font-bold text-blue-700"
                            )

                        unique_variants = len(
                            set(row.get("Variant", "") for row in validations_data)
                        )
                        with ui.column().classes("gap-1"):
                            ui.label("Unique Variants").classes("text-sm text-gray-600")
                            ui.label(str(unique_variants)).classes(
                                "text-3xl font-bold text-green-700"
                            )

                        unique_families = len(
                            set(row.get("FID", "") for row in validations_data)
                        )
                        with ui.column().classes("gap-1"):
                            ui.label("Families").classes("text-sm text-gray-600")
                            ui.label(str(unique_families)).classes(
                                "text-3xl font-bold text-purple-700"
                            )

                        unique_samples = len(
                            set(row.get("Sample", "") for row in validations_data)
                        )
                        with ui.column().classes("gap-1"):
                            ui.label("Samples").classes("text-sm text-gray-600")
                            ui.label(str(unique_samples)).classes(
                                "text-3xl font-bold text-orange-700"
                            )

                        # Show ignored count if any
                        if ignored_count > 0:
                            with ui.column().classes("gap-1"):
                                ui.label("Ignored").classes("text-sm text-gray-600")
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
                            row.get("Validation", "Unknown") for row in validations_data
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

                    # Inheritance Pattern Chart
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
                                },
                                "yAxis": {"type": "value", "name": "Count"},
                                "series": [
                                    {
                                        "type": "bar",
                                        "data": list(inheritance_counts.values()),
                                        "itemStyle": {"color": "#3b82f6"},
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
                            row.get("User", "Unknown") for row in validations_data
                        )
                        ui.echart(
                            {
                                "tooltip": {},
                                "xAxis": {
                                    "type": "category",
                                    "data": list(user_counts.keys()),
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

                    # Validations by Family
                    with ui.card().classes("flex-1 min-w-[400px] p-4"):
                        ui.label("Validations by Family").classes(
                            "text-lg font-semibold mb-2 text-gray-800"
                        )
                        family_counts = Counter(
                            row.get("FID", "Unknown") for row in validations_data
                        )
                        # Show top 10 families
                        top_families = dict(family_counts.most_common(10))
                        ui.echart(
                            {
                                "tooltip": {},
                                "yAxis": {
                                    "type": "category",
                                    "data": list(top_families.keys()),
                                },
                                "xAxis": {"type": "value", "name": "Count"},
                                "series": [
                                    {
                                        "type": "bar",
                                        "data": list(top_families.values()),
                                        "itemStyle": {"color": "#ec4899"},
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

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
