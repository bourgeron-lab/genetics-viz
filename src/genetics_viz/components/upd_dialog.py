"""Review dialog of one uniparental disomy: evidence, validation, diagnostic.

Opened from the family ROH tab (``upd_candidate_chrom``) and the family Sample
QC tab (a ``UPD_CANDIDATE`` flag). Both open the same record, keyed by family,
child and chromosome (:mod:`genetics_viz.utils.upd`).
"""

from pathlib import Path
from typing import Any, Callable, Dict, Optional

from nicegui import ui

from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.utils.auth import can_write, get_current_user
from genetics_viz.utils.diagnostic_files import diagnostic_path
from genetics_viz.utils.sample_qc import flag_title
from genetics_viz.utils.str_eh import diagnostic_badge
from genetics_viz.utils.upd import (
    ORIGINS,
    VALIDATIONS,
    default_origin,
    gather_evidence,
    latest_validation,
    normalise_chrom,
    read_upd_validations,
    save_upd_validation,
    upd_key,
)

#: Colour of a validation status chip.
VALIDATION_COLORS: Dict[str, str] = {
    "present": "#2e7d32",
    "absent": "#c62828",
    "uncertain": "#ef6c00",
}


def upd_status(
    data_dir: Path, family_id: str, sample: str, chrom: Any
) -> Dict[str, str]:
    """``{"validation": …, "origin": …, "diagnostic": …}`` of one UPD, ``""`` if none."""
    key = upd_key(chrom)
    latest = latest_validation(read_upd_validations(data_dir, family_id, key), sample)
    row: Dict[str, Any] = {}
    add_diagnostic_status_to_row(
        row,
        load_diagnostic_map(diagnostic_path(data_dir, "UPD"), family_id),
        key,
        sample,
    )
    return {
        "validation": (latest or {}).get("Validation", ""),
        "origin": (latest or {}).get("Inheritance", ""),
        "diagnostic": row.get("Diagnostic") or "",
    }


def status_chips(status: Dict[str, str]) -> None:
    """The validation and diagnostic of a UPD as two small chips."""
    validation = status.get("validation")
    if validation:
        origin = status.get("origin")
        ui.label(
            validation + (f" · {origin}" if origin and origin != "unknown" else "")
        ).classes("text-xs text-white rounded px-1").style(
            f"background-color: {VALIDATION_COLORS.get(validation, '#757575')}"
        ).tooltip("UPD validation")
    else:
        ui.label("not validated").classes("text-xs text-gray-400 italic")
    diagnostic = status.get("diagnostic")
    if diagnostic:
        badge = diagnostic_badge(diagnostic)
        ui.label(diagnostic).classes("text-xs text-white rounded px-1").style(
            f"background-color: {badge['color']}"
        ).tooltip("Curator diagnostic")


def show_upd_dialog(
    anchor: ui.element,
    data_dir: Path,
    family_id: str,
    sample: str,
    chrom: Any,
    *,
    on_saved: Optional[Callable[[], None]] = None,
) -> None:
    """Open the review of ``sample``'s ``chrom``; ``on_saved`` runs on close after a save.

    ``anchor`` must sit outside the caller's refreshable section: the save
    redraws that section, which would otherwise delete the open dialog.
    """
    chrom = normalise_chrom(chrom)
    key = upd_key(chrom)
    show_detail = can_write()
    evidence = gather_evidence(data_dir, family_id, sample, chrom)
    changed = {"value": False}

    def on_close(e: Any) -> None:
        if not e.value and changed["value"] and on_saved is not None:
            changed["value"] = False
            on_saved()

    anchor.clear()  # the previous review's dialog
    with anchor:
        dialog = ui.dialog()
    with dialog, ui.card().classes("w-[52rem] max-w-full"):
        dialog.on_value_change(on_close)
        ui.label(f"Uniparental disomy — {sample}, chr{chrom}").classes(
            "text-lg font-semibold text-blue-700"
        )

        ui.label("Evidence").classes("text-base font-semibold text-blue-700 mt-2")
        if evidence.roh:
            ui.label(
                f"ROH: an isolated block of long runs of homozygosity on chr{chrom},"
                " with little elsewhere: suggests uniparental isodisomy rather than"
                " parental relatedness."
            ).classes("text-sm")
        for f in evidence.flags:
            ui.label(
                f"Sample QC: {flag_title(f, show_detail)} ({' – '.join(f.samples)})"
                + (f": {f.evidence}" if show_detail and f.evidence else "")
            ).classes("text-sm").mark("upd-qc-evidence")
        if not evidence.roh and not evidence.flags:
            ui.label(
                "Neither the ROH nor the Sample QC tables point at this chromosome."
            ).classes("text-sm text-gray-500 italic")

        _render_validation(data_dir, family_id, sample, chrom, evidence, changed)
        _render_diagnostic(data_dir, family_id, sample, key, evidence.source, changed)

        with ui.row().classes("w-full justify-end mt-2"):
            ui.button("Close", on_click=dialog.close).props("flat")
    dialog.open()


def _render_validation(
    data_dir: Path,
    family_id: str,
    sample: str,
    chrom: str,
    evidence: Any,
    changed: Dict[str, bool],
) -> None:
    key = upd_key(chrom)
    ui.label("Validation").classes("text-base font-semibold text-blue-700 mt-2")

    @ui.refreshable
    def history() -> None:
        rows = [
            r
            for r in read_upd_validations(data_dir, family_id, key)
            if r.get("Sample") == sample
        ]
        if not rows:
            ui.label("No validation recorded for this chromosome.").classes(
                "text-sm text-gray-500 italic"
            )
            return
        ui.table(
            columns=[
                {"name": field, "label": label, "field": field, "align": "left"}
                for field, label in (
                    ("Validation", "Validation"),
                    ("Inheritance", "Origin"),
                    ("User", "User"),
                    ("Comment", "Comment"),
                    ("Date", "Date"),
                )
            ],
            rows=[
                {
                    **r,
                    "Date": (r.get("Timestamp") or "")[:10],
                    "Validation": r.get("Validation", "")
                    + (" (ignored)" if r.get("Ignore") == "1" else ""),
                    "id": i,
                }
                for i, r in enumerate(rows)
            ],
            row_key="id",
        ).props("dense flat").classes("w-full")

    history()

    if not can_write():
        ui.label("Curators and administrators can record validations.").classes(
            "text-xs text-gray-500"
        )
        return

    with ui.row().classes("items-end gap-3 w-full flex-wrap mt-2"):
        validation = (
            ui.select(VALIDATIONS, value="present", label="Validation")
            .props("outlined dense")
            .classes("w-36")
            .mark("upd-validation")
        )
        origin = (
            ui.select(
                ORIGINS,
                value=default_origin(evidence.flags, evidence.father, evidence.mother),
                label="Parental origin",
            )
            .props("outlined dense")
            .classes("w-40")
            .mark("upd-origin")
        )
        comment = ui.input("Comment").props("outlined dense").classes("flex-1 min-w-48")

        def save() -> None:
            try:
                save_upd_validation(
                    data_dir,
                    family_id=family_id,
                    sample=sample,
                    chrom=chrom,
                    user=get_current_user() or "",
                    validation=validation.value,
                    origin=origin.value,
                    comment=(comment.value or "").strip(),
                )
            except (PermissionError, ValueError) as ex:
                ui.notify(str(ex), type="negative")
                return
            except OSError as ex:
                ui.notify(f"Could not save the validation: {ex}", type="negative")
                return
            ui.notify(
                f"UPD validation saved: {sample} chr{chrom} {validation.value}",
                type="positive",
            )
            comment.value = ""
            changed["value"] = True
            history.refresh()

        ui.button("Save validation", icon="save", on_click=save).props(
            "color=primary"
        ).mark("upd-validation-save")
    ui.label(
        "Present or absent acknowledges the Sample QC flag (as confirmed or"
        " explained); uncertain leaves it counting."
    ).classes("text-xs text-gray-500")


def _render_diagnostic(
    data_dir: Path,
    family_id: str,
    sample: str,
    key: str,
    source: str,
    changed: Dict[str, bool],
) -> None:
    path = diagnostic_path(data_dir, "UPD")
    ui.label("Diagnostic").classes("text-base font-semibold text-blue-700 mt-2")
    # The review dialog lives here, outside the refreshable section: a refresh
    # from an ignore toggle would otherwise delete the dialog that sent it.
    anchor = ui.element("div")

    def open_review() -> None:
        with anchor:
            show_diagnostic_dialog(
                family_id=family_id,
                variant_key=key,
                gene="",
                impact=source,
                sample=sample,
                variant_type="upd",
                diagnostic_file=path,
                on_save_callback=on_change,
            )

    @ui.refreshable
    def section() -> None:
        row: Dict[str, Any] = {}
        add_diagnostic_status_to_row(
            row, load_diagnostic_map(path, family_id), key, sample
        )
        with ui.row().classes("items-center gap-2"):
            if row.get("Diagnostic"):
                badge = diagnostic_badge(row["Diagnostic"])
                ui.label(row["Diagnostic"]).classes(
                    "text-xs text-white rounded px-1"
                ).style(f"background-color: {badge['color']}")
            else:
                ui.label("none").classes("text-xs text-gray-400 italic")
            ui.button("Diagnostic review", icon="biotech", on_click=open_review).props(
                "flat dense no-caps"
            ).mark("upd-diagnostic-open")

    def on_change(_action: str) -> None:
        changed["value"] = True
        section.refresh()

    section()
