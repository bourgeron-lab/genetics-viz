"""Diagnostic review dialog for recording variant diagnostic conclusions."""

import csv
import getpass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from nicegui import ui

from genetics_viz.components.diagnostic_loader import (
    save_diagnostic_entry,
    update_diagnostic_ignore_status,
)


def show_diagnostic_dialog(
    family_id: str,
    variant_key: str,
    gene: str,
    impact: str,
    sample: str,
    variant_type: str,
    diagnostic_file: Path,
    on_save_callback: Optional[Callable[[str], None]] = None,
) -> None:
    """Show a lightweight diagnostic review dialog.

    Args:
        family_id: Family ID
        variant_key: Variant identifier (chr:pos:ref:alt or chr:start-end)
        gene: Gene symbol(s)
        impact: VEP consequence or SV call type (GAIN/LOSS)
        sample: Sample ID
        variant_type: "snv" or "sv"
        diagnostic_file: Path to the diagnostics TSV file
        on_save_callback: Optional callback invoked after save or ignore toggle
    """
    with (
        ui.dialog().props("full-width") as dialog,
        ui.card().classes("w-full max-w-3xl"),
    ):
        with ui.column().classes("w-full p-4 gap-4"):
            # ---- Header ----
            with ui.row().classes("items-center justify-between w-full"):
                ui.label("🔬 Diagnostic Review").classes(
                    "text-xl font-bold text-blue-900"
                )
                ui.button(icon="close", on_click=lambda: dialog.close()).props(
                    "flat round"
                )

            # ---- Variant summary ----
            with ui.card().classes("w-full bg-gray-50"):
                with ui.row().classes("gap-6 flex-wrap"):
                    _summary_field("Variant", variant_key)
                    _summary_field("Gene", gene)
                    _summary_field("Impact", impact)
                    _summary_field("Sample", sample)
                    _summary_field("Family", family_id)
                    _summary_field("Type", variant_type.upper())

            # ---- Form ----
            with ui.card().classes("w-full"):
                ui.label("New diagnostic").classes("font-semibold text-sm mb-2")
                with ui.row().classes("items-end gap-4 w-full flex-wrap"):
                    user_input = (
                        ui.input("User").props("outlined dense").classes("w-40")
                    )
                    user_input.value = getpass.getuser()

                    diagnostic_select = (
                        ui.select(
                            ["pathogenic", "uncertain", "benign"],
                            label="Diagnostic",
                        )
                        .props("outlined dense")
                        .classes("w-40")
                    )
                    diagnostic_select.value = "uncertain"

                    comment_input = (
                        ui.input("Comment")
                        .props("outlined dense")
                        .classes("flex-1 min-w-[200px]")
                    )

                    ui.button(
                        "Save",
                        icon="save",
                        on_click=lambda: _save_diagnostic(),
                    ).props("color=primary dense")

            # ---- History ----
            history_container = ui.column().classes("w-full")

            def load_history() -> None:
                """Load and display diagnostic history for this variant."""
                entries: list[Dict[str, Any]] = []

                if diagnostic_file.exists():
                    with open(diagnostic_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            if (
                                row.get("FID") == family_id
                                and row.get("Variant") == variant_key
                                and row.get("Sample") == sample
                            ):
                                entries.append(row)

                history_container.clear()

                with history_container:
                    ui.label("History").classes("font-semibold text-sm mb-1")
                    if not entries:
                        ui.label("No diagnostics recorded yet").classes(
                            "text-gray-500 text-sm italic"
                        )
                    else:
                        with ui.column().classes("gap-1 w-full"):
                            for entry in entries:
                                diag = entry.get("Diagnostic", "")
                                user = entry.get("User", "")
                                ts = entry.get("Timestamp", "")
                                comment = entry.get("Comment", "")
                                is_ignored = entry.get("Ignore", "0") == "1"

                                color, icon = _diag_style(diag)

                                row_cls = "items-center gap-2 w-full px-2 py-1"
                                if is_ignored:
                                    row_cls += " opacity-50"

                                with ui.row().classes(row_cls):
                                    ui.icon(icon).classes(f"text-{color} text-lg")
                                    ui.label(diag).classes(
                                        f"font-semibold text-{color} text-sm"
                                    )
                                    ui.label(f"by {user}").classes(
                                        "text-gray-600 text-xs"
                                    )
                                    if comment:
                                        ui.label(f"— {comment}").classes(
                                            "text-gray-500 text-xs italic"
                                        )
                                    ui.space()
                                    ui.label(_format_timestamp(ts)).classes(
                                        "text-gray-400 text-xs"
                                    )

                                    def make_ignore_handler(
                                        entry_ts: str,
                                    ) -> Callable:
                                        def handler(e: Any) -> None:
                                            new_val = "1" if e.value else "0"
                                            ok = update_diagnostic_ignore_status(
                                                diagnostic_file,
                                                family_id,
                                                variant_key,
                                                sample,
                                                entry_ts,
                                                new_val,
                                            )
                                            if ok:
                                                action = (
                                                    "ignored" if e.value else "restored"
                                                )
                                                ui.notify(
                                                    f"Diagnostic {action}",
                                                    type="info",
                                                )
                                                load_history()
                                                if on_save_callback:
                                                    on_save_callback(action)
                                            else:
                                                ui.notify(
                                                    "Failed to update",
                                                    type="negative",
                                                )

                                        return handler

                                    sw = ui.switch(
                                        "Ignore",
                                        value=is_ignored,
                                        on_change=make_ignore_handler(ts),
                                    ).classes("text-xs")
                                    if is_ignored:
                                        sw.props("color=grey")

            def _save_diagnostic() -> None:
                """Save the diagnostic entry."""
                user = user_input.value.strip()
                diag_value = diagnostic_select.value
                comment = comment_input.value.strip() if comment_input.value else ""

                if not user:
                    ui.notify("Please enter a username", type="warning")
                    return
                if not diag_value:
                    ui.notify("Please select a diagnostic", type="warning")
                    return

                timestamp = datetime.now().isoformat()

                try:
                    save_diagnostic_entry(
                        diagnostic_file=diagnostic_file,
                        family_id=family_id,
                        variant_key=variant_key,
                        gene=gene,
                        impact=impact,
                        sample=sample,
                        user=user,
                        timestamp=timestamp,
                        comment=comment,
                        diagnostic=diag_value,
                    )
                    ui.notify(f"Diagnostic saved: {diag_value}", type="positive")
                    dialog.close()
                    if on_save_callback:
                        on_save_callback(diag_value)
                except Exception as e:
                    ui.notify(f"Error saving diagnostic: {e}", type="negative")

            # Initial load
            load_history()

    dialog.open()


def _summary_field(label: str, value: str) -> None:
    """Render a label:value pair in the summary card."""
    with ui.column().classes("gap-0"):
        ui.label(label).classes("text-xs text-gray-500 uppercase")
        ui.label(value).classes("text-sm font-medium")


def _diag_style(diagnostic: str) -> tuple[str, str]:
    """Return (tailwind_color, material_icon) for a diagnostic value."""
    if diagnostic == "pathogenic":
        return "red", "error"
    elif diagnostic == "benign":
        return "green", "check_circle"
    elif diagnostic == "uncertain":
        return "orange", "help"
    elif diagnostic == "conflicting":
        return "amber", "warning"
    return "grey", "help_outline"


def _format_timestamp(ts: str) -> str:
    """Format an ISO timestamp to a shorter display form."""
    try:
        dt = datetime.fromisoformat(ts)
        return dt.strftime("%Y-%m-%d %H:%M")
    except (ValueError, TypeError):
        return ts
