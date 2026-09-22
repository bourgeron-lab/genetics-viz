"""The cohort's ghfc-ngs workflow parameters file, shown verbatim."""

from __future__ import annotations

from typing import Any, List

from nicegui import ui

from genetics_viz.utils.cohort_state import step_label
from genetics_viz.utils.pipeline_params import (
    get_requested_steps,
    load_params,
    read_params_text,
)


def render_params_panel(cohort: Any) -> None:
    """Render the cohort's parameters file.

    Shown verbatim rather than re-serialised from the parsed mapping: these
    files carry a great deal of their meaning in comments — which reference
    bundle a panel name pins, why a depth threshold was lowered — and none of
    that survives a round trip through the parser.
    """
    params_file = getattr(cohort, "params_file", None)
    if params_file is None:
        ui.label("No workflow parameters file for this cohort").classes(
            "text-gray-500 italic"
        )
        return

    text = read_params_text(params_file)
    steps: List[str] = get_requested_steps(load_params(params_file))

    with ui.card().classes("w-full"):
        # ui.code brings its own copy button, so the header carries only the
        # provenance: which file on disk this is.
        with ui.column().classes("gap-0 min-w-0 w-full"):
            ui.label("Workflow parameters").classes("text-lg font-semibold")
            # w-full because ui.column sets align-items: flex-start, so a
            # child is sized by its content and a long path escapes the card.
            # break-all rather than truncate: the path is provenance, and a
            # path has no spaces to wrap at.
            ui.label(str(params_file)).classes(
                "text-xs text-gray-400 font-mono break-all w-full"
            )

        if steps:
            with ui.row().classes("items-center gap-1 flex-wrap"):
                ui.label("Requested steps:").classes("text-xs text-gray-500")
                for step in steps:
                    ui.badge(step_label(step), color="blue").props("outline").classes(
                        "text-[10px]"
                    )

        if not text:
            ui.label("The file could not be read.").classes(
                "text-sm text-red-500 italic"
            )
            return

        with ui.element("div").classes(
            "w-full max-h-[70vh] overflow-auto border border-gray-200 rounded"
        ):
            ui.code(text, language="yaml").classes("w-full text-xs")
