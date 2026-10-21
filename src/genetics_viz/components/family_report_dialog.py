"""The family report: a dialog to read it and write the conclusion, and printing.

Both family pages put two buttons at the end of the title row, for curators
only. "See report" opens :func:`show_family_report`, "PDF" calls
:func:`print_family_report`. Either way the report is one self-contained HTML
document (:mod:`utils.report_html`), built off-thread:

- the dialog shows it in an ``<iframe srcdoc>``, so the page's Quasar and
  Tailwind styles stay out and the screen looks like the print;
- printing writes it into a hidden iframe and calls the browser's print. The
  document title is the default PDF file name.
"""

import asyncio
import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional

from nicegui import ui

from genetics_viz.components.offthread_panel import OffThreadPanel
from genetics_viz.utils.auth import can_write, get_current_user
from genetics_viz.utils.family_report import FamilyReport, build_family_report
from genetics_viz.utils.report_html import render_report_html
from genetics_viz.utils.report_summaries import (
    latest_family_summary,
    save_family_summary,
)

logger = logging.getLogger(__name__)


def _build(
    data_dir: Path, family_id: str, family: Any, cohort: Any, user: str
) -> Dict[str, Any]:
    report: FamilyReport = build_family_report(
        data_dir, family_id, family, cohort, user
    )
    return {"html": render_report_html(report), "title": report.document_title}


def print_html(html: str, title: str) -> None:
    """Print ``html`` through a hidden iframe; ``title`` names the PDF."""
    ui.run_javascript(
        f"""
(() => {{
  const frame = document.createElement('iframe');
  frame.style.cssText = 'position:fixed;right:0;bottom:0;width:0;height:0;border:0';
  frame.setAttribute('data-family-report-print', '');
  const pageTitle = document.title;
  frame.onload = () => {{
    const w = frame.contentWindow;
    document.title = {json.dumps(title)};
    w.addEventListener('afterprint', () => {{
      document.title = pageTitle;
      setTimeout(() => frame.remove(), 1000);
    }});
    w.focus();
    w.print();
  }};
  frame.srcdoc = {json.dumps(html)};
  document.body.appendChild(frame);
}})();
"""
    )


async def print_family_report(
    data_dir: Path, family_id: str, family: Any, cohort: Any = None
) -> None:
    """Build the report off-thread, then open the browser's print dialog."""
    if not can_write():
        return
    user = get_current_user()
    note = ui.notification("Preparing the report…", spinner=True, timeout=None)
    try:
        built = await asyncio.to_thread(
            _build, data_dir, family_id, family, cohort, user
        )
    except Exception as e:  # a bad file must not break the page
        logger.warning("Family report of %s failed", family_id, exc_info=True)
        note.dismiss()
        ui.notify(f"Could not build the report: {e}", type="negative")
        return
    note.dismiss()
    print_html(built["html"], built["title"])


def show_family_report(
    anchor: ui.element,
    data_dir: Path,
    family_id: str,
    family: Any,
    cohort: Any = None,
) -> None:
    """The report in a maximized dialog, with the curator conclusion editor."""
    if not can_write():
        return
    user = get_current_user()
    state: Dict[str, Optional[Dict[str, Any]]] = {"built": None}

    anchor.clear()  # the previous report's dialog
    with anchor:
        dialog = ui.dialog().props("maximized")
    with (
        dialog,
        ui.card().classes("w-full h-full p-0 gap-0").mark("family-report-dialog"),
    ):
        with ui.row().classes("w-full items-center gap-2 px-4 py-2 border-b"):
            ui.icon("description", color="primary").classes("text-2xl")
            ui.label(f"Report — family {family_id}").classes(
                "text-lg font-semibold text-blue-900"
            )
            ui.space()

            def do_print() -> None:
                if state["built"] is not None:
                    print_html(state["built"]["html"], state["built"]["title"])

            ui.button("Print / PDF", icon="picture_as_pdf", on_click=do_print).props(
                "unelevated color=primary no-caps"
            ).mark("family-report-print")
            ui.button(icon="close", on_click=dialog.close).props("flat round")

        with (
            ui.expansion("Curator conclusion", icon="edit_note")
            .classes("w-full px-2")
            .mark("family-report-conclusion")
        ):
            current = latest_family_summary(data_dir, family_id)
            text = (
                ui.textarea(
                    placeholder="What was found and whether recontacting the family is"
                    " advised. Printed at the top of the report.",
                    value=current.text if current else "",
                )
                .props("outlined autogrow")
                .classes("w-full")
                .mark("family-report-summary")
            )
            with ui.row().classes("w-full items-center gap-2"):
                saved_label = ui.label(
                    f"Last saved by {current.user}, {current.date}" if current else ""
                ).classes("text-xs text-gray-500")
                ui.space()

                def save() -> None:
                    try:
                        save_family_summary(
                            data_dir, family_id, text.value or "", user=user
                        )
                    except (PermissionError, ValueError, OSError) as e:
                        ui.notify(f"Not saved: {e}", type="negative")
                        return
                    saved = latest_family_summary(data_dir, family_id)
                    saved_label.text = (
                        f"Last saved by {saved.user}, {saved.date}" if saved else ""
                    )
                    ui.notify("Conclusion saved", type="positive")
                    reload()

                ui.button("Save conclusion", icon="save", on_click=save).props(
                    "unelevated no-caps"
                ).mark("family-report-save")

        panel = OffThreadPanel("family-report-body")
        panel.container.classes("flex-1 min-h-0")

        def render(built: Dict[str, Any]) -> None:
            state["built"] = built
            frame = ui.element("iframe").classes("w-full h-full border-0")
            frame.props["srcdoc"] = built["html"]
            frame.props["title"] = built["title"]
            frame.mark("family-report-frame")

        def reload() -> None:
            panel.load(
                lambda: _build(data_dir, family_id, family, cohort, user),
                render,
                "Building the report…",
            )

        reload()
    dialog.open()


def render_report_buttons(
    data_dir: Path, family_id: str, family: Any, cohort: Any = None
) -> None:
    """ "See report" and "PDF", for curators; nothing for viewers.

    The dialog hangs from an anchor created here, in the title row, outside
    every refreshable section of the page.
    """
    if not can_write():
        return
    anchor = ui.element()
    ui.button(
        "See report",
        icon="description",
        on_click=lambda: show_family_report(
            anchor, data_dir, family_id, family, cohort
        ),
    ).props("outline no-caps color=primary").mark("family-report-view")

    async def do_pdf() -> None:
        await print_family_report(data_dir, family_id, family, cohort)

    ui.button("PDF", icon="picture_as_pdf", on_click=do_pdf).props(
        "outline no-caps color=primary"
    ).mark("family-report-pdf").tooltip("Print the report or save it as PDF")
