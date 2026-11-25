"""Diagnostic review dialog for recording variant diagnostic conclusions."""

import csv
from genetics_viz.components.cell_line_note import render_cell_line_note
from genetics_viz.utils.auth import can_write, get_current_user
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, Mapping, Optional

from nicegui import ui

from genetics_viz.components.acmg_checklist import AcmgChecklist
from genetics_viz.components.diagnostic_loader import (
    save_diagnostic_entry,
    update_diagnostic_ignore_status,
)
from genetics_viz.utils import acmg
from genetics_viz.utils.diagnostic_badges import (
    UNCERTAIN,
    VERDICTS,
    verdict_icon,
    verdict_label,
)
from genetics_viz.utils.diagnostic_snapshot import SNAPSHOT_COLUMNS

#: Kinds whose dialog offers the ACMG criteria checklist: sequence variants,
#: and SVs for intragenic CNVs (the PVS1 tree). Whole-gene CNVs have their own
#: points (Riggs 2020), not covered yet; STR, UPD and aneuploidy have none.
ACMG_KINDS = ("snv", "sv")


def show_diagnostic_dialog(
    family_id: str,
    variant_key: str,
    gene: str,
    impact: str,
    sample: str,
    variant_type: str,
    diagnostic_file: Path,
    on_save_callback: Optional[Callable[[str], None]] = None,
    details: Optional[Mapping[str, str]] = None,
) -> None:
    """Show a lightweight diagnostic review dialog.

    Args:
        family_id: Family ID
        variant_key: Variant identifier (chr:pos:ref:alt, chr:start-end:type,
            or chr:pos-end:LOCUS for an STR)
        gene: Gene symbol(s)
        impact: VEP consequence, SV call type (GAIN/LOSS) or STR status
        sample: Sample ID
        variant_type: "snv", "sv", "str", "upd" or "aneuploidy"; "snv" and
            "sv" get the ACMG criteria checklist
        diagnostic_file: Path to the diagnostics TSV file
        on_save_callback: Optional callback invoked after save or ignore toggle
        details: The variant's snapshot from the table row
            (:mod:`genetics_viz.utils.diagnostic_snapshot`), saved with the
            diagnostic and shown under the summary. ``None`` saves none.
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

            # A cell line's de novo calls, UPDs and CNVs may come from culture.
            render_cell_line_note(family_id, [sample], "diagnostic")

            # ---- Variant summary ----
            with ui.card().classes("w-full bg-gray-50"):
                with ui.row().classes("gap-6 flex-wrap"):
                    _summary_field("Variant", variant_key)
                    _summary_field("Gene", gene)
                    _summary_field("Impact", impact)
                    _summary_field("Sample", sample)
                    _summary_field("Family", family_id)
                    _summary_field("Type", variant_type.upper())
                shown = [
                    (c, details[c])
                    for c in SNAPSHOT_COLUMNS
                    if c != "Source" and details and details.get(c)
                ]
                if shown:
                    ui.label(
                        f"Recorded with the diagnostic (from {details.get('Source')})"
                    ).classes("text-xs text-gray-500 uppercase mt-2")
                    with (
                        ui.row()
                        .classes("gap-x-6 gap-y-1 flex-wrap")
                        .mark("diagnostic-details")
                    ):
                        for column, value in shown:
                            _summary_field(column, value)

            # ---- Form ----
            with ui.card().classes("w-full"):
                ui.label("New diagnostic").classes("font-semibold text-sm mb-2")
                with ui.row().classes("items-end gap-4 w-full flex-wrap"):
                    user_input = (
                        ui.input("User")
                        .props("outlined dense readonly")
                        .classes("w-40")
                    )
                    user_input.value = get_current_user()

                    diagnostic_select = (
                        ui.select(
                            {v: verdict_label(v) for v in VERDICTS},
                            label="Diagnostic",
                            value=UNCERTAIN,
                            on_change=lambda: override_note.refresh(),
                        )
                        .props("outlined dense")
                        .classes("w-56")
                        .mark("diagnostic-verdict")
                    )

                    comment_input = (
                        ui.input("Comment")
                        .props("outlined dense")
                        .classes("flex-1 min-w-[200px]")
                    )

                    if can_write():
                        ui.button(
                            "Save",
                            icon="save",
                            on_click=lambda: _save_diagnostic(),
                        ).props("color=primary dense").mark("diagnostic-save")
                    else:
                        ui.label("Read-only access").classes(
                            "text-sm text-gray-400 italic"
                        )

                checklist: Optional[AcmgChecklist] = None
                if variant_type.lower() in ACMG_KINDS and can_write():

                    def follow_criteria(computed: str) -> None:
                        # The class the points give, until the curator picks
                        # another: then it is an override, said below.
                        if computed:
                            diagnostic_select.value = computed
                        override_note.refresh()

                    checklist = AcmgChecklist(on_change=follow_criteria)
                    checklist.render()

                @ui.refreshable
                def override_note() -> None:
                    if checklist is None or not checklist.applied:
                        return
                    computed = checklist.computed
                    if computed and computed != diagnostic_select.value:
                        ui.label(
                            f"The criteria give {verdict_label(computed).lower()}:"
                            f" saving {verdict_label(diagnostic_select.value).lower()}"
                            " records a curator override."
                        ).classes("text-xs text-amber-800").mark("acmg-override")

                override_note()

            # ---- History ----
            history_container = ui.column().classes("w-full")

            def load_history() -> None:
                """Load and display diagnostic history for this variant."""
                entries: list[Dict[str, Any]] = []

                if diagnostic_file.exists():
                    with open(diagnostic_file, "r") as f:
                        reader = csv.DictReader(
                            f, delimiter="\t", quoting=csv.QUOTE_NONE
                        )
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

                                icon, color = verdict_icon(diag)
                                icon = icon or "help_outline"

                                row_cls = "items-center gap-2 w-full px-2 py-1"
                                if is_ignored:
                                    row_cls += " opacity-50"

                                with ui.row().classes(row_cls):
                                    ui.icon(icon).classes(f"text-{color} text-lg")
                                    ui.label(verdict_label(diag)).classes(
                                        f"font-semibold text-{color} text-sm"
                                    )
                                    ui.label(f"by {user}").classes(
                                        "text-gray-600 text-xs"
                                    )
                                    if comment:
                                        ui.label(f"— {comment}").classes(
                                            "text-gray-500 text-xs italic"
                                        )
                                    if criteria := acmg.summary(entry, verdict_label):
                                        ui.label(criteria).classes(
                                            "text-xs text-gray-600 font-mono"
                                        ).mark("diagnostic-history-criteria")
                                    ui.space()
                                    ui.label(_format_timestamp(ts)).classes(
                                        "text-gray-400 text-xs"
                                    )

                                    def make_ignore_handler(
                                        entry_ts: str,
                                    ) -> Callable:
                                        def handler(e: Any) -> None:
                                            new_val = "1" if e.value else "0"
                                            try:
                                                ok = update_diagnostic_ignore_status(
                                                    diagnostic_file,
                                                    family_id,
                                                    variant_key,
                                                    sample,
                                                    entry_ts,
                                                    new_val,
                                                )
                                            except PermissionError as ex:
                                                ui.notify(str(ex), type="negative")
                                                return
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
                                    if not can_write():
                                        sw.props("disable")

            def _save_diagnostic() -> None:
                """Save the diagnostic entry."""
                if not can_write():
                    ui.notify("Permission denied", type="negative")
                    return
                user = get_current_user()
                diag_value = diagnostic_select.value
                comment = comment_input.value.strip() if comment_input.value else ""

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
                        details=details,
                        criteria=checklist.applied if checklist else (),
                    )
                    ui.notify(
                        f"Diagnostic saved: {verdict_label(diag_value)}",
                        type="positive",
                    )
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


def _format_timestamp(ts: str) -> str:
    """Format an ISO timestamp to a shorter display form."""
    try:
        dt = datetime.fromisoformat(ts)
        return dt.strftime("%Y-%m-%d %H:%M")
    except (ValueError, TypeError):
        return ts
