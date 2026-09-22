"""The ghfc-ngs parameters reference, shown beside a cohort's own params.yml."""

from __future__ import annotations

import json

from nicegui import ui

from genetics_viz.utils.ghfc_ngs_docs import (
    MARKDOWN_EXTRAS,
    PARAMS_DOC_PAGE,
    FetchedDoc,
    clear_doc_cache,
    fetch_params_doc,
)

#: ``min-w-0`` matters as much as the scrolling here: the document holds 18
#: tables, and without it a wide one's min-content floor propagates up the flex
#: chain and drags the whole row past the page. See the flex-sizing notes in
#: CLAUDE.md.
_SCROLL = "w-full min-w-0 max-h-[70vh] overflow-auto overflow-x-auto"

#: Marks the scroll container so the anchor handler below can find it from a
#: clicked link without depending on the DOM shape around it.
_SCROLL_MARKER = "data-ghfc-doc-scroll"

#: The document carries 73 in-page anchors. Left to the browser, clicking one
#: does scroll the panel to the heading -- but it also scrolls the whole page
#: (measured: window.scrollY 0 -> 995) and appends the fragment to the URL.
#: So: swallow the click, and move only the panel's own scrollTop by the
#: difference between the heading and the container. scrollIntoView would not
#: do, because it scrolls every scrollable ancestor too.
_ANCHOR_HANDLER = (
    """(e) => {
    const link = e.target.closest('a[href^="#"]');
    if (!link) return;
    e.preventDefault();
    const root = link.closest('[%s]');
    if (!root) return;
    const id = decodeURIComponent(link.getAttribute('href').slice(1));
    const target = Array.from(root.querySelectorAll('[id]')).find(n => n.id === id);
    if (!target) return;
    root.scrollTop +=
        target.getBoundingClientRect().top - root.getBoundingClientRect().top;
}"""
    % _SCROLL_MARKER
)


#: The document was written for a full-width GitHub page, so its default
#: heading sizes are absurd in a ~440px column — the h1 alone wrapped over three
#: lines. This scales the prose down to reference-material size and tightens the
#: tables, scoped to this panel so no other markdown in the app is affected.
_PROSE_CSS = """
[%(marker)s] .nicegui-markdown h1 { font-size: 1.25rem; margin: 0.6em 0 0.4em; }
[%(marker)s] .nicegui-markdown h2 { font-size: 1.05rem; margin: 1.1em 0 0.4em;
    border-bottom: 1px solid #e5e7eb; padding-bottom: 0.2em; }
[%(marker)s] .nicegui-markdown h3 { font-size: 0.925rem; margin: 1em 0 0.3em; }
[%(marker)s] .nicegui-markdown { font-size: 0.8125rem; line-height: 1.5; }
[%(marker)s] .nicegui-markdown table { font-size: 0.75rem; }
[%(marker)s] .nicegui-markdown th,
[%(marker)s] .nicegui-markdown td { padding: 0.2rem 0.4rem; }
[%(marker)s] .nicegui-markdown pre { font-size: 0.75rem; }
[%(marker)s] .nicegui-markdown code { font-size: 0.78rem; }
[%(marker)s] .nicegui-markdown blockquote { border-left: 3px solid #d1d5db;
    padding-left: 0.6rem; margin: 0.6em 0; color: #4b5563; }
""" % {"marker": _SCROLL_MARKER}


def _inject_prose_css_once() -> None:
    """Add the panel's prose styles once per client.

    Injected as a <style> element through run_javascript rather than with
    ui.add_css: this panel is built by a lazy tab's timer callback, long after
    the page head has been sent, so a head injection is simply dropped —
    measured, the rule never reached the document. tanstack_table.py notes the
    same trap for its own scripts.
    """
    client = ui.context.client
    if getattr(client, "_ghfc_doc_css_injected", False):
        return
    client._ghfc_doc_css_injected = True

    ui.run_javascript(f"""
        (function() {{
            if (document.getElementById('ghfc-doc-prose')) return;
            const style = document.createElement('style');
            style.id = 'ghfc-doc-prose';
            style.textContent = {json.dumps(_PROSE_CSS)};
            document.head.appendChild(style);
        }})();
    """)


def _github_link(label: str = "View on GitHub") -> None:
    """Link out to the rendered original."""
    with ui.link(target=PARAMS_DOC_PAGE, new_tab=True).classes(
        "text-xs text-blue-600 hover:underline no-underline whitespace-nowrap"
    ):
        with ui.row().classes("items-center gap-0.5 no-wrap"):
            ui.label(label)
            ui.icon("open_in_new", size="12px")


def render_params_doc_panel() -> None:
    """Render the upstream reference for the parameters file.

    Reads from the cache and stays synchronous, like the sibling panels; the
    caller warms the cache off-thread first, because the fetch is a network
    read and NiceGUI's event loop is single-threaded.
    """

    _inject_prose_css_once()

    @ui.refreshable
    def content() -> None:
        doc: FetchedDoc = fetch_params_doc()

        with ui.card().classes("w-full min-w-0"):
            with ui.row().classes(
                "w-full min-w-0 items-start justify-between no-wrap gap-2"
            ):
                with ui.column().classes("gap-0 min-w-0 w-full"):
                    ui.label("Parameters reference").classes("text-lg font-semibold")
                    with ui.row().classes("items-center gap-2 no-wrap min-w-0"):
                        ui.label(
                            "ghfc-ngs · fetched "
                            f"{doc.fetched_at.strftime('%Y-%m-%d %H:%M')}"
                        ).classes("text-xs text-gray-400 truncate min-w-0")
                        _github_link()

                def refresh() -> None:
                    clear_doc_cache()
                    content.refresh()

                ui.button(icon="refresh", on_click=refresh).props(
                    "flat dense round color=grey"
                ).tooltip("Re-fetch the document from GitHub")

            if not doc.ok:
                # No document is a degraded state, not an error: the rest of the
                # tab still works, and the original is one click away.
                with ui.row().classes("items-start gap-1 no-wrap w-full"):
                    ui.icon("cloud_off", size="14px").classes(
                        "text-amber-700 mt-0.5 shrink-0"
                    )
                    with ui.column().classes("gap-1 min-w-0 flex-1"):
                        ui.label(
                            "The parameters reference could not be loaded."
                        ).classes("text-xs text-amber-700")
                        ui.label(doc.error or "").classes(
                            "text-xs text-gray-500 break-all"
                        )
                        _github_link("Read it on GitHub instead")
                return

            # Delegated, so it works for markdown that arrives after the
            # container is built and survives a refresh.
            with (
                ui.element("div")
                .classes(_SCROLL)
                .props(_SCROLL_MARKER)
                .on("click", js_handler=_ANCHOR_HANDLER)
            ):
                ui.markdown(doc.markdown or "", extras=MARKDOWN_EXTRAS).classes(
                    "w-full min-w-0"
                )

    content()
