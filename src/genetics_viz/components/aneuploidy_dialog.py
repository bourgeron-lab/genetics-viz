"""Aneuploidy review: the SVs tab's card and the dialog of one sample's chromosome.

The card lists the family's whole-chromosome events (WisecondorX) and atypical
sex chromosomes (Sample QC karyotype), one row per member and chromosome
(:func:`genetics_viz.utils.aneuploidy.family_candidates`). Each row opens the
dialog: the evidence, the validation (present, absent or uncertain, with
direction, mosaicism, the ISCN karyotype and the origin) and the diagnostic.
The Sample QC tab opens the same dialog from a ``KARYOTYPE_<k>`` flag.

Readers see the sex chromosomes as "Sex chromosome anomaly", the Sample QC
reader wording, never the karyotype.
"""

from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from nicegui import ui

from genetics_viz.components.cell_line_note import render_cell_line_note
from genetics_viz.components.diagnostic_dialog import show_diagnostic_dialog
from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.utils import aneuploidy as an
from genetics_viz.utils.auth import can_write, get_current_user
from genetics_viz.utils.diagnostic_files import diagnostic_path
from genetics_viz.utils import sample_qc
from genetics_viz.utils.sample_qc import flag_title
from genetics_viz.utils.str_eh import diagnostic_badge

#: Colour of a validation status chip.
VALIDATION_COLORS: Dict[str, str] = {
    "present": "#2e7d32",
    "absent": "#c62828",
    "uncertain": "#ef6c00",
}


def _reader_sex_label() -> str:
    return str(sample_qc.KARYOTYPE_FLAG.get("reader_label") or "Sex chromosome anomaly")


def shown_label(candidate: an.AneuploidyCandidate, show_detail: bool) -> str:
    """``+21``, ``XXY``; for a reader, the sex chromosomes without the karyotype."""
    if candidate.chrom == an.SEX_CHROMS and not show_detail:
        return _reader_sex_label()
    return candidate.label


def shown_sources(candidate: an.AneuploidyCandidate, show_detail: bool) -> str:
    if show_detail or candidate.chrom != an.SEX_CHROMS:
        return " · ".join(candidate.sources)
    parts = [
        "WisecondorX " + ", ".join(e.text for e in candidate.events)
        if candidate.events
        else "",
        "Sample QC" if candidate.karyotype or candidate.flags else "",
    ]
    return " · ".join(p for p in parts if p)


def _evidence_source(candidate: an.AneuploidyCandidate) -> str:
    """The diagnostic's ``Impact``: which results point at it."""
    parts = []
    if candidate.events:
        parts.append("WisecondorX")
    if candidate.karyotype or candidate.flags:
        parts.append("Sample QC")
    return " + ".join(parts) or "Curator"


def _diagnostic(data_dir: Path, family_id: str, key: str, sample: str) -> str:
    row: Dict[str, Any] = {}
    add_diagnostic_status_to_row(
        row,
        load_diagnostic_map(diagnostic_path(data_dir, an.KIND), family_id),
        key,
        sample,
    )
    return row.get("Diagnostic") or ""


def status_chips(validation: str, diagnostic: str) -> None:
    """The validation and diagnostic as two small chips."""
    if validation:
        ui.label(validation).classes("text-xs text-white rounded px-1").style(
            f"background-color: {VALIDATION_COLORS.get(validation, '#757575')}"
        ).tooltip("Aneuploidy validation")
    else:
        ui.label("not validated").classes("text-xs text-gray-400 italic")
    if diagnostic:
        badge = diagnostic_badge(diagnostic)
        ui.label(diagnostic).classes("text-xs text-white rounded px-1").style(
            f"background-color: {badge['color']}"
        ).tooltip("Curator diagnostic")


def _mb(value: Any) -> str:
    try:
        return f"{float(value) / 1e6:,.1f}"
    except (TypeError, ValueError):
        return "—"


def _num(value: Any, digits: int = 2) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "—"


# ---------------------------------------------------------------------------
# The dialog
# ---------------------------------------------------------------------------


def show_aneuploidy_dialog(
    anchor: ui.element,
    data_dir: Path,
    family_id: str,
    sample: str,
    chrom: Any,
    members: Sequence[Tuple[str, str]],
    *,
    on_saved: Optional[Callable[[], None]] = None,
) -> None:
    """Open the review of ``sample``'s ``chrom`` (``"XY"`` for the sex chromosomes).

    ``members`` are the family's ``(sample, sex code)``
    (:func:`genetics_viz.utils.aneuploidy.member_sexes`). ``on_saved`` runs on
    close after a save. ``anchor`` must sit outside the caller's refreshable
    section: the save redraws that section, which would otherwise delete the
    open dialog.
    """
    chrom = an.record_chrom(chrom)
    show_detail = can_write()
    review = an.load_family_review(data_dir, family_id, members)
    declared = dict(members)
    candidate = review.candidate(sample, chrom) or an.AneuploidyCandidate(
        sample, chrom, declared=an.declared_karyotype(declared.get(sample, ""))
    )
    changed = {"value": False}

    def on_close(e: Any) -> None:
        if not e.value and changed["value"] and on_saved is not None:
            changed["value"] = False
            on_saved()

    anchor.clear()  # the previous review's dialog
    with anchor:
        dialog = ui.dialog()
    with dialog, ui.card().classes("w-[60rem] max-w-full").mark("aneuploidy-dialog"):
        dialog.on_value_change(on_close)
        ui.label(
            f"Aneuploidy — {sample}, {an.chrom_label(chrom)}"
            + (
                f" ({shown_label(candidate, show_detail)})"
                if shown_label(candidate, show_detail) != an.chrom_label(chrom)
                else ""
            )
        ).classes("text-lg font-semibold text-blue-700")
        # A cell line's CNVs may come from culture.
        render_cell_line_note(family_id, [sample], "cnv")

        ui.label("Evidence").classes("text-base font-semibold text-blue-700 mt-2")
        _render_evidence(review, candidate, show_detail)
        _render_validation(data_dir, family_id, candidate, changed, show_detail)
        _render_diagnostic(data_dir, family_id, candidate, changed)

        with ui.row().classes("w-full justify-end mt-2"):
            ui.button("Close", on_click=dialog.close).props("flat").mark(
                "aneuploidy-dialog-close"
            )
    dialog.open()


def _render_evidence(
    review: an.FamilyReview, candidate: an.AneuploidyCandidate, show_detail: bool
) -> None:
    sample, chrom = candidate.sample, candidate.chrom
    if candidate.events:
        ui.label(
            "WisecondorX: "
            + ", ".join(e.text for e in candidate.events)
            + " of the chromosome's assessable span (the q arm of an acrocentric"
            " one) in calls past the permissive thresholds"
        ).classes("text-sm").mark("aneuploidy-wcx-evidence")
    segments = review.segments(sample, chrom)
    if segments:
        ui.table(
            columns=[
                {"name": k, "label": label, "field": k, "align": "left"}
                for k, label in (
                    ("chr", "Chromosome"),
                    ("start", "Start (Mb)"),
                    ("end", "End (Mb)"),
                    ("ratio", "log2 ratio"),
                    ("zscore", "z-score"),
                    ("type", "Call"),
                )
            ],
            rows=[
                {
                    "id": i,
                    "chr": r.get("chr", ""),
                    "start": _mb(r.get("start")),
                    "end": _mb(r.get("end")),
                    "ratio": _num(r.get("ratio")),
                    "zscore": _num(r.get("zscore"), 1),
                    # The classification the WisecondorX table shows.
                    "type": r.get("wisecondorX") or r.get("source") or "",
                }
                for i, r in enumerate(segments)
            ],
            row_key="id",
        ).props("dense flat").classes("w-full").mark("aneuploidy-segments")
    qc = review.qc_rows.get(sample)
    if chrom == an.SEX_CHROMS and qc is not None:
        if show_detail:
            declared = f"; declared {candidate.declared}" if candidate.declared else ""
            ui.label(
                f"Sample QC: karyotype {qc.get('karyotype') or 'unknown'} (X"
                f" {_num(qc.get('x_copies'))}, Y {_num(qc.get('y_copies'))}"
                f" copies){declared}"
            ).classes("text-sm").mark("aneuploidy-qc-evidence")
        elif candidate.karyotype:
            ui.label(f"Sample QC: {_reader_sex_label().lower()}").classes("text-sm")
    for f in candidate.flags:
        ui.label(
            f"Sample QC flag: {flag_title(f, show_detail)}"
            + (f": {f.evidence}" if show_detail and f.evidence else "")
        ).classes("text-sm").mark("aneuploidy-flag-evidence")
    if not candidate.events and not segments and not candidate.karyotype:
        if not candidate.flags:
            ui.label(
                "Neither WisecondorX nor Sample QC point at this chromosome."
            ).classes("text-sm text-gray-500 italic")
    ui.label(
        "Each call opens in the WisecondorX table's review, with its copy-number"
        " view and read evidence."
    ).classes("text-xs text-gray-500")


def _render_validation(
    data_dir: Path,
    family_id: str,
    candidate: an.AneuploidyCandidate,
    changed: Dict[str, bool],
    show_detail: bool,
) -> None:
    sample, chrom, key = candidate.sample, candidate.chrom, candidate.key
    hide_karyotype = chrom == an.SEX_CHROMS and not show_detail
    ui.label("Validation").classes("text-base font-semibold text-blue-700 mt-2")

    @ui.refreshable
    def history() -> None:
        rows = [
            r
            for r in an.read_aneuploidy_validations(data_dir, family_id, key)
            if r.get("Sample") == sample
        ]
        if not rows:
            ui.label("No validation recorded for this chromosome.").classes(
                "text-sm text-gray-500 italic"
            )
            return
        fields = [
            ("Validation", "Validation"),
            ("Direction", "Direction"),
            ("Mosaic", "Mosaic"),
            ("Karyotype", "Karyotype"),
            ("Inheritance", "Origin"),
            ("User", "User"),
            ("Comment", "Comment"),
            ("Date", "Date"),
        ]
        if hide_karyotype:
            fields = [f for f in fields if f[0] != "Karyotype"]
        ui.table(
            columns=[
                {"name": f, "label": label, "field": f, "align": "left"}
                for f, label in fields
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
        ).props("dense flat").classes("w-full").mark("aneuploidy-validation-history")

    history()

    if not show_detail:
        ui.label("Curators and administrators can record validations.").classes(
            "text-xs text-gray-500"
        )
        return

    with ui.row().classes("items-end gap-3 w-full flex-wrap mt-2"):
        validation = (
            ui.select(an.VALIDATIONS, value="present", label="Validation")
            .props("outlined dense")
            .classes("w-36")
            .mark("aneuploidy-validation")
        )
        direction = (
            ui.select(
                {d: d or "—" for d in an.DIRECTIONS},
                value=candidate.direction
                if candidate.direction in an.DIRECTIONS
                else "",
                label="Direction",
            )
            .props("outlined dense")
            .classes("w-32")
            .mark("aneuploidy-direction")
        )
        mosaic = (
            ui.select({m: m or "—" for m in an.MOSAIC}, value="", label="Mosaic")
            .props("outlined dense")
            .classes("w-28")
            .mark("aneuploidy-mosaic")
        )
        karyotype = (
            ui.input("Karyotype (ISCN)", value=candidate.suggested_karyotype)
            .props("outlined dense")
            .classes("w-44")
            .mark("aneuploidy-karyotype")
        )
        origin = (
            ui.select(an.ORIGINS, value="unknown", label="Origin")
            .props("outlined dense")
            .classes("w-36")
            .mark("aneuploidy-origin")
        )
    with ui.row().classes("items-end gap-3 w-full flex-wrap"):
        comment = ui.input("Comment").props("outlined dense").classes("flex-1 min-w-48")

        def save() -> None:
            try:
                an.save_aneuploidy_validation(
                    data_dir,
                    family_id=family_id,
                    sample=sample,
                    chrom=chrom,
                    user=get_current_user() or "",
                    validation=validation.value,
                    origin=origin.value,
                    direction=direction.value or "",
                    mosaic=mosaic.value or "",
                    karyotype=(karyotype.value or "").strip(),
                    comment=(comment.value or "").strip(),
                )
            except (PermissionError, ValueError) as ex:
                ui.notify(str(ex), type="negative")
                return
            except OSError as ex:
                ui.notify(f"Could not save the validation: {ex}", type="negative")
                return
            ui.notify(
                f"Aneuploidy validation saved: {sample} {an.chrom_label(chrom)}"
                f" {validation.value}",
                type="positive",
            )
            comment.value = ""
            changed["value"] = True
            history.refresh()

        ui.button("Save validation", icon="save", on_click=save).props(
            "color=primary"
        ).mark("aneuploidy-validation-save")
    if chrom == an.SEX_CHROMS:
        ui.label(
            "Present or absent acknowledges the sample's Sample QC karyotype flag"
            " (as confirmed or explained); uncertain leaves it counting."
        ).classes("text-xs text-gray-500")


def _render_diagnostic(
    data_dir: Path,
    family_id: str,
    candidate: an.AneuploidyCandidate,
    changed: Dict[str, bool],
) -> None:
    path = diagnostic_path(data_dir, an.KIND)
    sample, key = candidate.sample, candidate.key
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
                impact=_evidence_source(candidate),
                sample=sample,
                variant_type="aneuploidy",
                diagnostic_file=path,
                on_save_callback=on_change,
            )

    @ui.refreshable
    def section() -> None:
        diagnostic = _diagnostic(data_dir, family_id, key, sample)
        with ui.row().classes("items-center gap-2"):
            if diagnostic:
                badge = diagnostic_badge(diagnostic)
                ui.label(diagnostic).classes("text-xs text-white rounded px-1").style(
                    f"background-color: {badge['color']}"
                )
            else:
                ui.label("none").classes("text-xs text-gray-400 italic")
            ui.button("Diagnostic review", icon="biotech", on_click=open_review).props(
                "flat dense no-caps"
            ).mark("aneuploidy-diagnostic-open")

    def on_change(_action: str) -> None:
        changed["value"] = True
        section.refresh()

    section()


# ---------------------------------------------------------------------------
# The SVs tab's card
# ---------------------------------------------------------------------------


def render_aneuploidy_card(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
    on_saved: Optional[Callable[[], None]] = None,
    listeners: Optional[List[Callable[[], None]]] = None,
) -> None:
    """The family's aneuploidy candidates, for the selected members.

    Draws nothing when the family has none. ``on_saved`` runs after a review
    saved something (the page's header badge and member icons follow a
    karyotype flag's acknowledgement). ``listeners`` are run after a review
    anywhere on the page (the Sample QC tab opens the same dialog); the card
    adds its own reload.
    """
    data_dir = store.data_dir
    members = an.member_sexes(store, family_id, cohort_name)
    show_detail = can_write()
    state = {"review": an.load_family_review(data_dir, family_id, members)}
    # Outside the refreshable section, which a save redraws.
    anchor = ui.element("div")

    def reload() -> None:
        state["review"] = an.load_family_review(data_dir, family_id, members)
        card.refresh()

    def saved() -> None:
        reload()
        if on_saved is not None:
            on_saved()

    def open_review(sample: str, chrom: str) -> None:
        show_aneuploidy_dialog(
            anchor, data_dir, family_id, sample, chrom, members, on_saved=saved
        )

    @ui.refreshable
    def card() -> None:
        selected = set(selected_members["value"])
        shown = [c for c in state["review"].candidates if c.sample in selected]
        if not shown:
            return
        with ui.card().classes("w-full mb-2").mark("aneuploidy-card"):
            with ui.row().classes("items-center gap-2"):
                ui.label("Aneuploidy").classes("text-base font-semibold text-blue-700")
                to_review = an.count_to_review(shown)
                if show_detail and to_review:
                    ui.label(f"{to_review} to review").classes(
                        "text-xs text-white rounded px-1 bg-orange-600"
                    ).mark("aneuploidy-to-review")
                ui.label(
                    "Whole-chromosome gains and losses (WisecondorX) and atypical"
                    " sex chromosomes (Sample QC)"
                ).classes("text-xs text-gray-500")
            for c in shown:
                with (
                    ui.row()
                    .classes("items-center gap-3 w-full border-b border-gray-100 py-1")
                    .mark(f"aneuploidy-row-{c.sample}-{c.chrom}")
                ):
                    ui.label(c.sample).classes("text-sm font-medium w-24")
                    ui.label(shown_label(c, show_detail)).classes(
                        "text-sm font-semibold w-40"
                    )
                    ui.label(shown_sources(c, show_detail)).classes(
                        "text-xs text-gray-600 flex-1 min-w-0 truncate"
                    )
                    status_chips(
                        c.status, _diagnostic(data_dir, family_id, c.key, c.sample)
                    )

                    def make_handler(sample: str, chrom: str) -> Callable[[], None]:
                        return lambda: open_review(sample, chrom)

                    ui.button(
                        "Review",
                        icon="rate_review",
                        on_click=make_handler(c.sample, c.chrom),
                    ).props("flat dense no-caps").mark(
                        f"aneuploidy-review-{c.sample}-{c.chrom}"
                    )

    card()
    data_table_refreshers.append(card.refresh)
    if listeners is not None:
        listeners.append(reload)
