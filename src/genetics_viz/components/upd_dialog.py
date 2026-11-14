"""Review dialog of one uniparental disomy: evidence, chromosome view, validation, diagnostic.

Opened from the family UPD tab (one row per child and chromosome), the ROH tab
(``upd_candidate_chrom``) and the Sample QC tab (a ``UPD_CANDIDATE`` flag).
All open the same record, keyed by family, child and chromosome
(:mod:`genetics_viz.utils.upd`).

The chromosome view draws the upd step's ``windows.tsv`` (from ``upd_v4``):
the alleles only one parent has and the Mendelian errors against each, the
het rate and the depth, over lanes of the step's segments, the ROH runs and
the Sample QC span (:func:`genetics_viz.utils.upd_plots.chromosome_option`).
"""

from pathlib import Path
from typing import Any, Callable, Dict, Optional

from nicegui import ui

from genetics_viz.components.cell_line_note import render_cell_line_note
from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.offthread_panel import OffThreadPanel
from genetics_viz.utils import upd_calls
from genetics_viz.utils.auth import can_write, get_current_user
from genetics_viz.utils.diagnostic_files import diagnostic_path
from genetics_viz.utils.sample_qc import flag_title
from genetics_viz.utils.str_eh import diagnostic_badge
from genetics_viz.utils.upd import (
    COMPOSITIONS,
    EXTENTS,
    ORIGINS,
    VALIDATIONS,
    UpdEvidence,
    default_origin,
    gather_evidence,
    latest_validation,
    normalise_chrom,
    read_upd_validations,
    roh_runs,
    save_upd_validation,
    span_from_evidence,
    upd_key,
    validated_span,
)
from genetics_viz.utils.upd_plots import chromosome_height, chromosome_option

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


def _fmt(value: Any, digits: int = 2) -> str:
    try:
        return f"{float(value):,.{digits}f}"
    except (TypeError, ValueError):
        return "—"


def _render_events(evidence: UpdEvidence) -> None:
    """The upd step's events on this chromosome, every confidence."""
    if not evidence.events:
        return
    ui.label(f"UPD step ({evidence.upd_label})").classes("text-sm font-medium")
    with ui.column().classes("w-full gap-1").mark("upd-step-evidence"):
        for e in evidence.events:
            info = upd_calls.confidence_info(e.get("confidence"))
            origin = upd_calls.origin_info(e.get("origin"))
            with ui.row().classes("items-center gap-2 flex-wrap text-sm"):
                ui.label(info["label"]).classes(
                    "text-xs text-white rounded px-1"
                ).style(f"background-color: {info['color']}")
                ui.label(upd_calls.span_text(e.get("start"), e.get("end"))).classes(
                    "font-medium"
                )
                ui.label(origin["label"]).style(f"color: {origin['color']}")
                ui.label(upd_calls.event_type(e))
                mode = upd_calls.mode_info(e.get("mode"))
                ui.label(f"({mode['label']})").classes("text-gray-500").tooltip(
                    mode["help"]
                )
                ui.label(
                    f"LOD {_fmt(e.get('lod'), 0)} · {e.get('n_strong') or 0:,} strong"
                    f" sites · support {_fmt(e.get('support_fraction'))}"
                    f" · depth {_fmt(e.get('depth_ratio'))}"
                ).classes("text-xs text-gray-600")
                if e.get("recurrent") == "yes":
                    ui.label("recurrent").classes(
                        "text-xs text-white rounded px-1 bg-gray-500"
                    ).tooltip(
                        f"Overlaps events of {e.get('recurrent_families')} other"
                        " families: most likely an artefact region"
                    )
                if e.get("shared_with"):
                    ui.label(f"shared with {e['shared_with']}").classes(
                        "text-xs text-white rounded px-1 bg-gray-500"
                    ).tooltip(
                        "A sibling has an overlapping event: a parental CNV or"
                        " mismapped reads, not a UPD"
                    )
            if e.get("imprinting"):
                ui.label(f"Imprinting: {e['imprinting']}").classes(
                    "text-xs text-purple-800 ml-2"
                )
            if e.get("notes"):
                ui.label(e["notes"]).classes("text-xs text-gray-500 ml-2")


def _render_chromosome_view(
    data_dir: Path,
    family_id: str,
    sample: str,
    chrom: str,
    evidence: UpdEvidence,
) -> None:
    ui.label("Chromosome view").classes("text-base font-semibold text-blue-700 mt-2")
    panel = OffThreadPanel("upd-chromosome-view")
    qc_span = next(
        (s for s in (span_from_evidence(f) for f in evidence.flags) if s), None
    )
    validated = validated_span(
        latest_validation(
            read_upd_validations(data_dir, family_id, upd_key(chrom)), sample
        )
    )

    def load() -> Dict[str, Any]:
        calls = upd_calls.load_family_upd(data_dir, family_id, evidence.upd_label)
        windows = []
        segments = []
        if calls is not None:
            windows = upd_calls.windows_of(upd_calls.load_windows(calls), sample, chrom)
            segments = calls.segments_of(sample, chrom)
        return {
            "windows": windows,
            "segments": segments,
            "roh": roh_runs(data_dir, family_id, sample, chrom),
            "label": calls.label if calls is not None else "",
        }

    def render(result: Dict[str, Any]) -> None:
        if not result["windows"]:
            ui.label(
                "No per-window genotype counts: the upd step wrote them from"
                " upd_v4 on"
                + (
                    f" (these tables are {result['label']})."
                    if result["label"]
                    else "."
                )
                + " Only the segments, ROH runs and Sample QC span are drawn."
            ).classes("text-xs text-gray-500")
        else:
            ui.label(
                "Up: alleles only the mother has (father homozygous, child"
                " heterozygous) and, darker, sites where the child lacks the"
                " father's allele. Down: the same for the father. A maternal UPD"
                " has no paternal-only allele; an isodisomy has no het call; a"
                " deletion of one copy looks like an isodisomy at half the depth."
            ).classes("text-xs text-gray-500")
        option = chromosome_option(
            chrom,
            result["windows"],
            segments=result["segments"],
            roh_runs=result["roh"],
            qc_span=qc_span,
            validated=validated,
            events=evidence.events,
        )
        ui.echart(option).classes("w-full").style(
            f"height: {chromosome_height(bool(result['windows']))}px"
        ).mark("upd-chromosome-chart")

    panel.load(load, render, "Loading the chromosome view…")


def show_upd_dialog(
    anchor: ui.element,
    data_dir: Path,
    family_id: str,
    sample: str,
    chrom: Any,
    *,
    on_saved: Optional[Callable[[], None]] = None,
    on_closed: Optional[Callable[[], None]] = None,
    upd_label: Optional[str] = None,
) -> None:
    """Open the review of ``sample``'s ``chrom``; ``on_saved`` runs on close after a save.

    ``on_closed`` runs on every close, after ``on_saved``.
    ``anchor`` must sit outside the caller's refreshable section: the save
    redraws that section, which would otherwise delete the open dialog.
    ``upd_label`` is the upd step's label shown (default: the newest).
    """
    chrom = normalise_chrom(chrom)
    key = upd_key(chrom)
    show_detail = can_write()
    evidence = gather_evidence(data_dir, family_id, sample, chrom, upd_label)
    changed = {"value": False}

    def on_close(e: Any) -> None:
        if e.value:
            return
        if changed["value"] and on_saved is not None:
            changed["value"] = False
            on_saved()
        if on_closed is not None:
            on_closed()

    anchor.clear()  # the previous review's dialog
    with anchor:
        dialog = ui.dialog()
    with dialog, ui.card().classes("w-[72rem] max-w-full").mark("upd-dialog"):
        dialog.on_value_change(on_close)
        ui.label(f"Uniparental disomy — {sample}, chr{chrom}").classes(
            "text-lg font-semibold text-blue-700"
        )
        render_cell_line_note(family_id, [sample], "upd")

        ui.label("Evidence").classes("text-base font-semibold text-blue-700 mt-2")
        _render_events(evidence)
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
        if not evidence.roh and not evidence.flags and not evidence.events:
            ui.label(
                "Neither the upd step, the ROH nor the Sample QC tables point at"
                " this chromosome."
            ).classes("text-sm text-gray-500 italic")

        _render_chromosome_view(data_dir, family_id, sample, chrom, evidence)
        _render_validation(data_dir, family_id, sample, chrom, evidence, changed)
        _render_diagnostic(data_dir, family_id, sample, key, evidence.source, changed)

        with ui.row().classes("w-full justify-end mt-2"):
            ui.button("Close", on_click=dialog.close).props("flat").mark(
                "upd-dialog-close"
            )
    dialog.open()


def _type_text(row: Dict[str, str]) -> str:
    parts = [
        upd_calls.composition_label(row.get("Composition")),
        upd_calls.extent_label(row.get("Extent")),
    ]
    return " · ".join(p for p in parts if p)


def _span_text(row: Dict[str, str]) -> str:
    span = validated_span(row)
    return upd_calls.span_text(*span) if span else ""


def _render_validation(
    data_dir: Path,
    family_id: str,
    sample: str,
    chrom: str,
    evidence: UpdEvidence,
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
                    ("Type", "Type"),
                    ("Span", "Span"),
                    ("UpdLabel", "Checked on"),
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
                    "Type": _type_text(r),
                    "Span": _span_text(r),
                    "id": i,
                }
                for i, r in enumerate(rows)
            ],
            row_key="id",
        ).props("dense flat").classes("w-full").mark("upd-validation-history")

    history()

    if not can_write():
        ui.label("Curators and administrators can record validations.").classes(
            "text-xs text-gray-500"
        )
        return

    event = evidence.event or {}
    qc_span = next(
        (s for s in (span_from_evidence(f) for f in evidence.flags) if s), None
    )
    if event.get("start") is not None and event.get("end") is not None:
        span = (int(event["start"]), int(event["end"]))
    else:
        span = qc_span
    composition = (
        event.get("composition") if event.get("composition") in COMPOSITIONS else ""
    )
    extent = event.get("extent") if event.get("extent") in EXTENTS else ""

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
                value=default_origin(
                    evidence.flags, evidence.father, evidence.mother, evidence.events
                ),
                label="Parental origin",
            )
            .props("outlined dense")
            .classes("w-40")
            .mark("upd-origin")
        )
        composition_select = (
            ui.select(
                {c: upd_calls.composition_label(c) or "—" for c in COMPOSITIONS},
                value=composition,
                label="Composition",
            )
            .props("outlined dense")
            .classes("w-40")
            .mark("upd-composition")
        )
        extent_select = (
            ui.select(
                {c: upd_calls.extent_label(c) or "—" for c in EXTENTS},
                value=extent,
                label="Extent",
            )
            .props("outlined dense")
            .classes("w-44")
            .mark("upd-extent")
        )
        start_mb = (
            ui.number(
                "Start (Mb)",
                value=round(span[0] / 1e6, 3) if span else None,
                format="%.3f",
                min=0,
                step=0.001,
            )
            .props("outlined dense")
            .classes("w-32")
            .mark("upd-start")
        )
        end_mb = (
            ui.number(
                "End (Mb)",
                value=round(span[1] / 1e6, 3) if span else None,
                format="%.3f",
                min=0,
                step=0.001,
            )
            .props("outlined dense")
            .classes("w-32")
            .mark("upd-end")
        )
    with ui.row().classes("items-end gap-3 w-full flex-wrap"):
        comment = ui.input("Comment").props("outlined dense").classes("flex-1 min-w-48")

        def save() -> None:
            start = start_mb.value
            end = end_mb.value
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
                    start=None if start is None else int(round(float(start) * 1e6)),
                    end=None if end is None else int(round(float(end) * 1e6)),
                    composition=composition_select.value or "",
                    extent=extent_select.value or "",
                    upd_label=evidence.upd_label,
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
        " explained) and settles the ROH candidate; uncertain leaves both"
        " counting. The span and type are what you validated: the report shows"
        " them even after a newer upd version calls the event differently."
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
