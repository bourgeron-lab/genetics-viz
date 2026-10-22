"""SV dialog: the validation form (curated position included) and its history."""

import csv
import traceback
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from nicegui import ui

from genetics_viz.components.validation_loader import (
    SV_VALIDATION_HEADER,
    parents_to_validate,
    validated_samples,
)
from genetics_viz.utils.auth import can_write, get_current_user
from genetics_viz.utils.tsv import append_tsv_rows, clean_tsv_field

_INHERITANCE_OPTIONS = [
    "unknown",
    "de novo",
    "paternal",
    "maternal",
    "not paternal",
    "not maternal",
    "either",
    "homozygous",
]

_VALIDATION_OPTIONS = ["present", "absent", "uncertain", "different"]

_HISTORY_BADGE_COLORS = {
    "present": "bg-green-600",
    "absent": "bg-red-600",
    "uncertain": "bg-yellow-600",
    "different": "bg-orange-600",
}


@dataclass
class SVValidationForm:
    """The inputs other panels fill: the curated position and inheritance."""

    curated_start: Any
    curated_end: Any
    inheritance: Any

    def curated_input(self, key: str) -> Any:
        return self.curated_start if key == "start" else self.curated_end


def render_sv_validation(
    *,
    validation_file: Path,
    family_id: str,
    variant_key: str,
    chrom: str,
    start: str,
    end: str,
    sample: str,
    sample_parents: Dict[str, Optional[str]],
    bounds: Dict[str, int],
    on_bounds_changed: Callable[[], None],
    on_saved: Callable[[], None],
) -> SVValidationForm:
    """The "SV Validation" form, then the sample's previous validations.

    The refresh button next to the curated position moves ``bounds`` to it
    and calls ``on_bounds_changed``. A save records the session user; with
    ``present`` it also records the parents the inheritance names, in one
    locked write. ``on_saved`` runs after a save.
    """
    ui.label("SV Validation").classes("text-xl font-semibold mb-2 mt-4")

    with ui.card().classes("w-full p-4 mb-4"):
        with ui.column().classes("gap-4"):
            with ui.row().classes("items-center gap-4 w-full flex-wrap"):
                ui.label("User:").classes("font-semibold")
                user_input = (
                    ui.input("Username")
                    .props("outlined dense readonly")
                    .classes("w-48")
                )
                user_input.value = get_current_user()

                ui.label("Inheritance:").classes("font-semibold ml-4")
                inheritance_select = (
                    ui.select(_INHERITANCE_OPTIONS, value="unknown")
                    .props("outlined dense")
                    .classes("w-40")
                )

                ui.label("Validation:").classes("font-semibold ml-4")
                validation_select = (
                    ui.select(_VALIDATION_OPTIONS, value="present")
                    .props("outlined dense")
                    .classes("w-40")
                )

                ui.label("Curated Position:").classes("font-semibold")
                curated_start_input = (
                    ui.input("Curated start").props("outlined dense").classes("w-48")
                )

                ui.label("-").classes("font-semibold ml-4")
                curated_end_input = (
                    ui.input("Curated end").props("outlined dense").classes("w-48")
                )

                def _refresh_roi_from_curated():
                    """Update ROI from curated position inputs."""
                    new_start = (
                        curated_start_input.value.strip()
                        if curated_start_input.value
                        else ""
                    )
                    new_end = (
                        curated_end_input.value.strip()
                        if curated_end_input.value
                        else ""
                    )
                    if not new_start and not new_end:
                        ui.notify(
                            "Enter curated start and/or end first",
                            type="warning",
                        )
                        return
                    try:
                        if new_start:
                            bounds["start"] = int(new_start)
                        if new_end:
                            bounds["end"] = int(new_end)
                    except ValueError:
                        ui.notify(
                            "Curated positions must be integers",
                            type="warning",
                        )
                        return
                    on_bounds_changed()
                    ui.notify(
                        f"ROI updated: {chrom}:{bounds['start']}-{bounds['end']}",
                        type="positive",
                    )

                ui.button(
                    icon="refresh",
                    on_click=_refresh_roi_from_curated,
                ).props("flat round dense size=sm color=blue").tooltip(
                    "Refresh IGV ROI to curated coordinates"
                )

            with ui.row().classes("items-center gap-4 w-full"):
                ui.label("Comment:").classes("font-semibold")
                comment_input = (
                    ui.input("Optional comment")
                    .props("outlined dense")
                    .classes("flex-grow")
                )

                if can_write():
                    ui.button(
                        "Save Validation",
                        icon="save",
                        on_click=lambda: save_sv_validation(),
                    ).props("color=blue")
                else:
                    ui.label("Read-only access").classes("text-sm text-gray-400 italic")

            def save_sv_validation():
                """Save an SV validation."""
                if not can_write():
                    ui.notify("Permission denied", type="negative")
                    return
                user = get_current_user()
                inheritance = inheritance_select.value
                validation_status = validation_select.value
                curated_start = (
                    curated_start_input.value.strip()
                    if curated_start_input.value
                    else ""
                )
                curated_end = (
                    curated_end_input.value.strip() if curated_end_input.value else ""
                )
                comment = clean_tsv_field(
                    comment_input.value.strip() if comment_input.value else ""
                )

                timestamp = datetime.now().isoformat()

                try:
                    # A validated present call names the parent(s) it came
                    # from: record them as carriers too, unless they are
                    # validated already. One locked write.
                    parent_samples_added = []
                    if validation_status == "present":
                        parent_samples_added = parents_to_validate(
                            inheritance,
                            sample_parents,
                            validated_samples(validation_file, family_id, variant_key),
                        )
                    curated = [curated_start, curated_end]
                    rows = [
                        [family_id, variant_key, sample, user, inheritance]
                        + [validation_status, *curated, comment, "0"]
                        + [timestamp]
                    ] + [
                        [family_id, variant_key, parent_id, user, "unknown"]
                        + ["present", *curated, f"(inherited from {sample})"]
                        + ["0", timestamp]
                        for parent_id in parent_samples_added
                    ]
                    append_tsv_rows(validation_file, SV_VALIDATION_HEADER, rows)

                    # Notify user about saved validations
                    notify_message = f"SV validation saved: {validation_status}"
                    if parent_samples_added:
                        notify_message += (
                            f" (also saved for {', '.join(parent_samples_added)})"
                        )
                    ui.notify(notify_message, type="positive")

                    # Reload validation history
                    load_sv_validation_history()

                    on_saved()

                except Exception as e:
                    ui.notify(f"Error saving validation: {e}", type="negative")
                    print(traceback.format_exc())

    # Validation history container
    validation_history_container = ui.column().classes("w-full mb-4")

    def load_sv_validation_history():
        """Load and display SV validation history."""
        validations = []

        if validation_file.exists():
            with open(validation_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
                for row in reader:
                    if (
                        row.get("FID") == family_id
                        and row.get("Variant") == variant_key
                        and row.get("Sample") == sample
                    ):
                        validations.append(row)

        validation_history_container.clear()

        with validation_history_container:
            ui.label("Previous validations:").classes("font-semibold mb-2")
            if not validations:
                ui.label("No validations recorded yet").classes(
                    "text-gray-500 text-sm italic"
                )
                return
            with ui.card().classes("w-full p-2"):
                with ui.column().classes("gap-2"):
                    for validation in validations:
                        _render_history_row(validation, start, end)

    load_sv_validation_history()

    return SVValidationForm(
        curated_start=curated_start_input,
        curated_end=curated_end_input,
        inheritance=inheritance_select,
    )


def _render_history_row(validation: Dict[str, str], start: str, end: str) -> None:
    val_status = validation.get("Validation", "")
    user = validation.get("User", "")
    timestamp = validation.get("Timestamp", "")
    comment = validation.get("Comment", "")
    curated_start = validation.get("CuratedStart", "")
    curated_end = validation.get("CuratedEnd", "")
    is_ignored = validation.get("Ignore", "0") == "1"

    # Format timestamp
    try:
        dt = datetime.fromisoformat(timestamp)
        formatted_time = dt.strftime("%Y-%m-%d %H:%M:%S")
    except Exception:
        formatted_time = timestamp

    # Color based on status
    if is_ignored:
        badge_color = "bg-gray-400"
    else:
        badge_color = _HISTORY_BADGE_COLORS.get(val_status, "bg-gray-500")

    with ui.row().classes("items-center gap-2 p-2 w-full"):
        ui.badge(val_status).classes(f"{badge_color} text-xs")
        ui.label(user).classes("text-sm font-semibold")
        ui.label(formatted_time).classes("text-xs text-gray-500")
        if curated_start or curated_end:
            ui.label(f"Curated: {curated_start or start}-{curated_end or end}").classes(
                "text-xs text-blue-600"
            )
        if comment:
            ui.label(f'"{comment}"').classes("text-sm text-gray-600 italic")
        if is_ignored:
            ui.badge("IGNORED").classes("bg-gray-400 text-xs")
