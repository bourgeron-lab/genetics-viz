"""Pipeline run status for a cohort, as shown on the home card and cohort page.

One module for both, so the two views cannot drift. The numbers come from the
cohort's ``.ghfc-ngs.state.json`` and are never computed here -- see
:mod:`genetics_viz.utils.cohort_state` for why.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from nicegui import ui

from genetics_viz.utils.cohort_state import (
    CohortState,
    RunRecord,
    StepRow,
    build_step_rows,
    clear_sha_cache,
    clear_state_cache,
    file_drift,
    format_duration,
    format_timestamp,
    read_cohort_state,
    record_status_meta,
    status_meta,
    summarize_steps,
)
from genetics_viz.utils.pipeline_params import (
    get_requested_steps,
    load_params,
)

#: Shown wherever the pipeline has left no record. Worth spelling out: the
#: schema documentation is explicit that absence is not evidence the cohort was
#: never processed.
_NO_RECORD_NOTE = (
    "No run recorded for this cohort. The pipeline writes "
    ".ghfc-ngs.state.json when it runs, so its absence does not mean the "
    "cohort has never been processed."
)

_BAR_TRACK = "h-1.5 rounded bg-gray-200 overflow-hidden w-full"


def _params_steps(cohort: Any) -> List[str]:
    """The steps the cohort's parameters file asks for."""
    params_file = getattr(cohort, "params_file", None)
    if params_file is None:
        return []
    return get_requested_steps(load_params(params_file))


def _pct_color(row: StepRow) -> str:
    """Bar colour by how far along a step is."""
    pct = row.completion.pct or 0.0
    if row.complete:
        return "#16a34a"  # green-600
    if pct >= 75:
        return "#2563eb"  # blue-600
    if pct >= 25:
        return "#d97706"  # amber-600
    return "#dc2626"  # red-600


def _render_step_row(row: StepRow, *, compact: bool = False) -> None:
    """One step: label, marker, bar and counts.

    Every text here can shrink: ``min-w-0`` on each flex item and ``truncate``
    rather than ``whitespace-nowrap``. A flex item defaults to
    ``min-width: auto`` and refuses to go below its content, and that
    min-content contribution propagates all the way up -- one unshrinkable
    label was enough to push the whole home-page card past its own width.
    """
    with ui.column().classes("w-full min-w-0 gap-0.5"):
        with ui.row().classes(
            "w-full min-w-0 items-baseline justify-between gap-2 no-wrap"
        ):
            with ui.row().classes("items-baseline gap-1.5 no-wrap min-w-0"):
                ui.label(row.label).classes("text-xs text-gray-700 truncate min-w-0")
                # A step the last run was not asked to do, or one the params
                # file has since dropped: the number, if any, does not describe
                # what the pipeline is currently being asked for. Dropped on the
                # compact card, which has no room for it.
                if not compact:
                    if row.in_params and not row.requested:
                        ui.label("not in last run").classes(
                            "text-[10px] text-gray-400 italic truncate min-w-0"
                        )
                    elif row.requested and not row.in_params:
                        ui.label("no longer requested").classes(
                            "text-[10px] text-gray-400 italic truncate min-w-0"
                        )
                if row.note:
                    with ui.icon("info", size="12px").classes("text-gray-300 shrink-0"):
                        ui.tooltip(row.note).classes("text-xs max-w-xs")
            if row.measured:
                ui.label(
                    f"{row.completion.pct:.0f}% · {row.detail(compact=compact)}"
                ).classes("text-[11px] text-gray-500 truncate min-w-0 text-right")
            else:
                # Never a 0% bar: unmeasured and zero are different facts.
                ui.label("not measured").classes(
                    "text-[11px] text-gray-400 italic truncate min-w-0 text-right"
                )

        with ui.element("div").classes(_BAR_TRACK):
            if row.measured:
                pct = max(0.0, min(100.0, row.completion.pct or 0.0))
                ui.element("div").classes("h-full rounded").style(
                    f"width: {pct:.1f}%; background-color: {_pct_color(row)}"
                )


def _render_step_rows(rows: List[StepRow], *, compact: bool = False) -> None:
    """The step table, or a note when there are no steps to show at all."""
    if not rows:
        ui.label("No pipeline steps declared in the parameters file.").classes(
            "text-xs text-gray-500 italic"
        )
        return
    with ui.column().classes("w-full min-w-0 gap-2"):
        for row in rows:
            _render_step_row(row, compact=compact)


def _note(text: str, *, tone: str = "gray", icon: Optional[str] = None) -> None:
    """A one-line caveat under the status header."""
    colors = {
        "gray": "text-gray-500",
        "amber": "text-amber-700",
        "red": "text-red-600",
        "green": "text-green-700",
    }
    with ui.row().classes("items-start gap-1 no-wrap w-full"):
        if icon:
            ui.icon(icon, size="14px").classes(f"{colors[tone]} mt-0.5 shrink-0")
        # min-w-0 flex-1: a flex item defaults to min-width:auto and will not
        # shrink below its content, so a long caveat would widen the panel.
        ui.label(text).classes(f"text-xs {colors[tone]} min-w-0 flex-1")


def _render_drift(cohort: Any, record: RunRecord) -> None:
    """Flag inputs edited since the run that produced these outputs.

    Hashes the local files: the paths inside the record are compute-cluster
    paths and do not resolve here.
    """
    pedigree = file_drift(
        record.pedigree_sha256, getattr(cohort, "pedigree_file", None)
    )
    params = file_drift(record.params_sha256, getattr(cohort, "params_file", None))

    if pedigree == "changed":
        _note(
            "The pedigree has been edited since this run — the counts below "
            "describe a different pedigree.",
            tone="amber",
            icon="warning",
        )
    if params == "changed":
        _note(
            "The parameters file has been edited since this run.",
            tone="amber",
            icon="warning",
        )
    if pedigree == "unchanged" and params in ("unchanged", "unknown"):
        _note("Pedigree unchanged since this run.", tone="green", icon="check")


def _render_run_header(state: CohortState) -> None:
    """Status chip, timings and every caveat that applies to the numbers."""
    record = state.last_run
    meta = status_meta(state)

    with ui.row().classes("items-center gap-2 flex-wrap"):
        ui.badge(meta.get("label", ""), color=meta.get("color", "grey")).classes(
            "text-xs"
        )
        if record is not None:
            ui.label(
                f"{format_timestamp(record.started_at)} → "
                f"{format_timestamp(record.finished_at)}"
            ).classes("text-xs text-gray-600")
            if record.duration is not None:
                ui.label(f"({format_duration(record.duration)})").classes(
                    "text-xs text-gray-400"
                )

    if record is None:
        return

    if state.last_successful_run is not None:
        ui.label(
            f"Last completed: {format_timestamp(state.last_successful_run.finished_at)}"
        ).classes("text-xs text-gray-600")
    else:
        ui.label("Last completed: never").classes("text-xs text-gray-500")

    if record.pipeline_version or record.short_commit:
        parts = [p for p in (record.pipeline_version, record.short_commit) if p]
        ui.label("Pipeline " + " @ ".join(parts)).classes("text-xs text-gray-400")

    if record.reason:
        _note(record.reason, tone="red", icon="error_outline")

    if record.is_stale_running:
        _note(
            "This run has reported nothing for over a day. There is no "
            "heartbeat, so a run killed by a walltime limit stays 'running' "
            "until the next run closes it out — it has most likely died.",
            tone="amber",
            icon="warning",
        )

    if record.measured_before_run:
        _note(
            "These counts were taken before the run started, so they do not "
            "include anything it produced.",
            tone="amber",
            icon="schedule",
        )

    if record.outputs_may_be_incomplete:
        _note(
            "This run failed. Nextflow cancels in-flight output copies when a "
            "run aborts, so the counts may under-report what was produced.",
            tone="amber",
            icon="warning",
        )


def _render_history(records: List[RunRecord]) -> None:
    """Previous terminal runs, newest first."""
    if not records:
        return
    with (
        ui.expansion(f"Run history ({len(records)})")
        .classes("w-full text-sm")
        .props("dense")
    ):
        with ui.grid(columns=4).classes("gap-x-3 gap-y-1 text-xs w-full"):
            hdr = "text-gray-500 font-semibold border-b border-gray-200 pb-1"
            for title in ("Status", "Started", "Duration", "Pipeline"):
                ui.label(title).classes(hdr)
            for rec in records:
                meta = record_status_meta(rec)
                ui.badge(
                    meta.get("label", rec.status), color=meta.get("color", "grey")
                ).classes("text-[10px] justify-self-start")
                ui.label(format_timestamp(rec.started_at)).classes("text-gray-600")
                ui.label(format_duration(rec.duration)).classes("text-gray-600")
                ui.label(
                    " @ ".join(p for p in (rec.pipeline_version, rec.short_commit) if p)
                    or "—"
                ).classes("text-gray-400 truncate")
                if rec.reason:
                    ui.label(rec.reason).classes("col-span-4 text-red-500 italic pb-1")


# --------------------------------------------------------------------------
# Public entry points
# --------------------------------------------------------------------------


def status_headline(cohort: Any) -> Dict[str, Any]:
    """Cheap summary for the collapsed home-card expansion.

    One stat plus one small JSON read, both cached, so this is safe to call for
    every card on the home page.
    """
    state = read_cohort_state(cohort.path)
    meta = status_meta(state)
    if state is not None and state.unsupported:
        return {
            **meta,
            "caption": f"state file schema {state.schema_version} not supported",
        }
    rows = build_step_rows(state, _params_steps(cohort))
    summary = summarize_steps(rows)
    return {**meta, "caption": summary.text, "all_complete": summary.all_complete}


def render_status_summary(cohort: Any) -> None:
    """Per-step rows only: the body of the home card's Status expansion."""
    state = read_cohort_state(cohort.path)
    if state is not None and state.unsupported:
        _note(
            f"This cohort's run record announces schema version "
            f"{state.schema_version}, which this build cannot read.",
            tone="amber",
            icon="warning",
        )
        return
    if state is None:
        _note(_NO_RECORD_NOTE)
    _render_step_rows(build_step_rows(state, _params_steps(cohort)), compact=True)


def render_status_panel(
    cohort: Any, *, title: Optional[str] = "Pipeline status"
) -> None:
    """Full run record for a cohort.

    ``title`` is None when the caller already names the panel — the Workflow
    tab puts it in an expansion header that carries the same title, chip and
    summary, and repeating it inside would read as two headings.
    """

    @ui.refreshable
    def content() -> None:
        state = read_cohort_state(cohort.path)

        with ui.card().classes("w-full"):
            with ui.row().classes("w-full items-center justify-between no-wrap"):
                if title:
                    ui.label(title).classes("text-lg font-semibold")
                else:
                    # Keeps the refresh button pushed to the right.
                    ui.element("div")

                def refresh() -> None:
                    clear_state_cache()
                    clear_sha_cache()
                    content.refresh()

                ui.button(icon="refresh", on_click=refresh).props(
                    "flat dense round color=grey"
                ).tooltip("Re-read the run record from disk")

            if state is not None and state.unsupported:
                _note(
                    f"This cohort's run record announces schema version "
                    f"{state.schema_version}; this build reads version 1. "
                    "Nothing from it is shown, because a field may have "
                    "changed meaning.",
                    tone="amber",
                    icon="warning",
                )
                ui.label(str(state.path)).classes(
                    "text-xs text-gray-400 font-mono break-all w-full"
                )
                return

            if state is None:
                _note(_NO_RECORD_NOTE)
            else:
                _render_run_header(state)
                if state.last_run is not None:
                    _render_drift(cohort, state.last_run)

            ui.separator()
            _render_step_rows(build_step_rows(state, _params_steps(cohort)))

            if state is not None and state.history:
                ui.separator()
                _render_history(state.history)

            if state is not None:
                ui.label(str(state.path)).classes(
                    "text-xs text-gray-400 font-mono break-all w-full mt-1"
                )

    content()
