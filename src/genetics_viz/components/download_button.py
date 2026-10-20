"""The Download button of the cohort and family tabs.

One export gives a plain button; several give a menu with one entry per TSV
plus, optionally, "All tables (.zip)". Each export is built when it is
clicked, from the tab's state at that moment, off the event loop: a builder
must not touch ``app.storage`` (resolve the user's role before, at render).

The content always reaches ``ui.download`` as bytes. NiceGUI's test client
takes a ``str`` for a URL and fetches it.
"""

import asyncio
import logging
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional

from nicegui import ui

from genetics_viz.utils.export import zip_texts

logger = logging.getLogger(__name__)

_TSV = "text/tab-separated-values"


@dataclass
class Download:
    """One TSV a tab offers."""

    #: Short id, used in the menu entry's test mark.
    key: str
    label: str
    filename: Callable[[], str]
    #: The TSV text; ``None``, ``""`` or a header alone when there is nothing
    #: to export.
    build: Callable[[], Optional[str]]


def _has_rows(text: Optional[str]) -> bool:
    return bool(text) and text.count("\n") > 1


async def _build(download: Download) -> Optional[str]:
    try:
        return await asyncio.to_thread(download.build)
    except Exception:
        logger.warning("Failed to build export %s", download.key, exc_info=True)
        ui.notify(f"Could not build {download.label}", type="negative")
        return None


async def _send(download: Download) -> None:
    text = await _build(download)
    if not _has_rows(text):
        ui.notify("Nothing to download", type="info")
        return
    ui.download(text.encode("utf-8"), download.filename(), media_type=_TSV)


async def _send_zip(downloads: List[Download], filename: Callable[[], str]) -> None:
    files: Dict[str, str] = {}
    for download in downloads:
        text = await _build(download)
        if _has_rows(text):
            files[download.filename()] = text
    if not files:
        ui.notify("Nothing to download", type="info")
        return
    ui.download(zip_texts(files), filename(), media_type="application/zip")


def download_button(
    downloads: List[Download],
    mark: str,
    zip_filename: Optional[Callable[[], str]] = None,
) -> ui.button:
    """``Download``, a button for one export or a menu for several.

    ``mark`` tags the button; menu entries are tagged ``<mark>-<key>`` and the
    zip entry ``<mark>-zip``.
    """
    props = "outline color=blue size=sm"
    if len(downloads) == 1:
        only = downloads[0]
        return (
            ui.button("Download", icon="download", on_click=lambda: _send(only))
            .props(props)
            .tooltip(only.label)
            .mark(mark)
        )

    def entry(download: Download) -> Callable:
        async def handler() -> None:
            await _send(download)

        return handler

    with ui.button("Download", icon="download").props(props).mark(mark) as button:
        with ui.menu():
            for download in downloads:
                ui.menu_item(download.label, on_click=entry(download)).mark(
                    f"{mark}-{download.key}"
                )
            if zip_filename is not None:
                ui.separator()
                ui.menu_item(
                    "All tables (.zip)",
                    on_click=lambda: _send_zip(downloads, zip_filename),
                ).mark(f"{mark}-zip")
    return button
