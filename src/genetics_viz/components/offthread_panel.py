"""A container whose content is computed off-thread and redrawn on demand.

Shared by the SV dialog's copy-number and read-evidence panels: show a
spinner, run a blocking loader in a thread, then render, unless a newer
refresh started meanwhile (adding a sample while the previous list loads).
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any, Callable

from nicegui import ui

logger = logging.getLogger(__name__)


class OffThreadPanel:
    def __init__(self, mark: str) -> None:
        self._generation = 0
        self.container = ui.column().classes("w-full gap-1").mark(mark)

    def show_message(self, text: str, classes: str = "text-gray-500 p-4") -> None:
        self._generation += 1
        self.container.clear()
        with self.container:
            ui.label(text).classes(classes)

    def load(
        self,
        loader: Callable[[], Any],
        render: Callable[[Any], None],
        loading_text: str,
    ) -> None:
        """Run ``loader`` in a thread, then ``render`` its result in the container."""
        self._generation += 1
        generation = self._generation
        self.container.clear()
        with self.container:
            with ui.row().classes("items-center gap-2 p-4"):
                ui.spinner(size="md")
                ui.label(loading_text).classes("text-gray-500")

            async def run() -> None:
                try:
                    result = await asyncio.to_thread(loader)
                except Exception as e:  # a bad file must not break the dialog
                    logger.warning("%s failed: %s", loading_text, e, exc_info=True)
                    if generation == self._generation:
                        self.show_message(f"Could not load: {e}", "text-red-600 p-4")
                    return
                if generation != self._generation:
                    return  # superseded
                self.container.clear()
                with self.container:
                    render(result)

            ui.timer(0.1, run, once=True)
