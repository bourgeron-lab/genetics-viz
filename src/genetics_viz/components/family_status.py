"""The family status on the family pages: a chip for everyone, a dialog for curators.

The chip shows the status and, when there is one, a reason to look at the
family again (``Reopened`` after a pedigree edit, ``Check`` otherwise), with
the details in its tooltip. It is drawn at once from the status file; the
reasons (pedigree edits, diagnostics, QC flags, step outputs) are read off
the event loop just after, and again after each save or curation change.
"""

import asyncio
import logging
from pathlib import Path
from typing import Any, Callable, List, Optional, Sequence

from nicegui import ui

from genetics_viz.utils import family_status as fs
from genetics_viz.utils.auth import can_write, get_current_user

logger = logging.getLogger(__name__)


class FamilyStatusChip:
    """The status chip of one family, in the page's title row."""

    def __init__(
        self,
        data_dir: Path,
        family_id: str,
        family: Any,
        cohorts: Sequence[Any] = (),
        on_saved: Optional[Callable[[], None]] = None,
    ) -> None:
        self.data_dir = data_dir
        self.family_id = family_id
        self.family = family
        self.cohorts = list(cohorts)
        self.on_saved = on_saved
        self.curator = can_write()
        #: The family's step outputs, walked once per page.
        self._results: Optional[List[Any]] = None
        self.state = fs.FamilyState(
            family_id,
            fs.family_statuses(data_dir).get(family_id),
            affected=fs.has_affected(family) if family is not None else True,
        )
        # The dialog hangs from here, outside every refreshable section.
        self._anchor = ui.element()
        with ui.row().classes("items-center gap-1 no-wrap") as self._row:

            @ui.refreshable
            def draw() -> None:
                self._render()

            draw()
        self._draw = draw
        ui.timer(0.2, self.update, once=True)

    # -- drawing -----------------------------------------------------------

    def _render(self) -> None:
        state = self.state
        info = fs.status_info(state.status)
        tooltip = fs.status_tooltip(state, self.curator)
        if self.curator:
            tooltip += "\nClick to set the status"
        chip = (
            ui.badge(f"Status: {info['label']}", color=info["color"])
            .classes("text-sm px-2 py-1")
            .tooltip(tooltip)
            .mark("family-status")
        )
        if self.curator:
            chip.classes("cursor-pointer").on("click", self.open_dialog)
        if state.reasons:
            reason = (
                ui.badge("Reopened" if state.reopened else "Check", color="#b45309")
                .classes("text-sm px-2 py-1")
                .tooltip("\n".join(r.text for r in state.reasons))
                .mark("family-status-reason")
            )
            if self.curator:
                reason.classes("cursor-pointer").on("click", self.open_dialog)

    async def update(self) -> None:
        """Read the family's state (with every reason) in a thread, then redraw."""
        try:
            if self._results is None:
                self._results = await asyncio.to_thread(
                    fs.step_outputs, self.data_dir, self.family_id
                )
            self.state = await asyncio.to_thread(
                fs.load_family_state,
                self.data_dir,
                self.family_id,
                self.family,
                self.cohorts,
                self._results,
            )
        except Exception:
            logger.warning("Failed to read %s's status", self.family_id, exc_info=True)
        self._draw.refresh()

    def refresh(self) -> None:
        """:meth:`update` from a sync callback (a save, a review closed).

        The timer hangs from the chip's row, which outlives any dialog.
        """
        with self._row:
            ui.timer(0.05, self.update, once=True)

    # -- dialog ------------------------------------------------------------

    def open_dialog(self) -> None:
        if not can_write():
            return
        state = self.state
        history = fs.family_history(self.data_dir, self.family_id)
        self._anchor.clear()
        with self._anchor:
            dialog = ui.dialog()
        with (
            dialog,
            ui.card().classes("w-[36rem] max-w-full").mark("family-status-dialog"),
        ):
            with ui.row().classes("w-full items-center"):
                ui.label(f"Status — family {self.family_id}").classes(
                    "text-lg font-semibold text-blue-900"
                )
                ui.space()
                ui.button(icon="close", on_click=dialog.close).props("flat round")

            record = state.record
            info = fs.status_info(state.status)
            with ui.row().classes("items-center gap-2"):
                ui.badge(info["label"], color=info["color"])
                ui.label(
                    f"{record.user}, {record.date}" if record else "No status set yet"
                ).classes("text-sm text-gray-600")
            if record and record.comment:
                ui.label(record.comment).classes("text-sm text-gray-700 italic")

            if state.reasons:
                with (
                    ui.column()
                    .classes(
                        "w-full gap-1 p-2 rounded bg-amber-50 border border-amber-300"
                    )
                    .mark("family-status-reasons")
                ):
                    for r in state.reasons:
                        with ui.row().classes("items-start gap-2 no-wrap"):
                            ui.icon("warning", color="amber-9").classes("text-base")
                            ui.label(r.text).classes("text-sm text-amber-900")

            options = {k: v["label"] for k, v in fs.STATUSES.items()}
            current = state.status if state.status in options else fs.NOT_REVIEWED
            select = (
                ui.select(options, value=current, label="Status")
                .props("outlined dense")
                .classes("w-full")
                .mark("family-status-select")
            )
            description = ui.label(
                fs.status_info(current).get("description", "")
            ).classes("text-xs text-gray-500")
            select.on_value_change(
                lambda e: description.set_text(
                    fs.status_info(e.value).get("description", "")
                )
            )
            comment = (
                ui.input(
                    label="Comment",
                    placeholder="What was found, or why (shown to curators only)",
                )
                .props("outlined dense")
                .classes("w-full")
                .mark("family-status-comment")
            )
            if state.waiting_rerun:
                ui.label(
                    "The pedigree was edited and ghfc-ngs has not run since:"
                    " the results still read the older pedigree."
                ).classes("text-xs text-amber-800")

            async def save() -> None:
                user = get_current_user()
                try:
                    await fs.save_family_status_async(
                        self.data_dir,
                        self.family_id,
                        select.value,
                        comment.value or "",
                        user=user or "",
                    )
                except (PermissionError, ValueError) as e:
                    ui.notify(str(e), type="negative")
                    return
                except OSError as e:
                    logger.warning("Failed to save a status", exc_info=True)
                    ui.notify(f"Could not save the status: {e}", type="negative")
                    return
                ui.notify(
                    f"{self.family_id}: {fs.status_label(select.value)}",
                    type="positive",
                )
                dialog.close()
                self.refresh()
                if self.on_saved is not None:
                    self.on_saved()

            with ui.row().classes("w-full justify-end gap-2"):
                ui.button("Cancel", on_click=dialog.close).props("flat no-caps")
                ui.button("Save", icon="save", on_click=save).props(
                    "unelevated no-caps color=primary"
                ).mark("family-status-save")

            if history:
                with (
                    ui.expansion(f"History ({len(history)})", icon="history")
                    .classes("w-full")
                    .mark("family-status-history")
                ):
                    for r in history:
                        with ui.row().classes("items-baseline gap-2 no-wrap"):
                            ui.label(r.date).classes(
                                "text-xs text-gray-500 whitespace-nowrap"
                            )
                            ui.badge(
                                fs.status_label(r.status),
                                color=fs.status_info(r.status)["color"],
                            ).classes("text-xs")
                            ui.label(r.user).classes("text-xs text-gray-600")
                            if r.comment:
                                ui.label(r.comment).classes("text-xs text-gray-700")
        dialog.open()


def render_family_status(
    data_dir: Path,
    family_id: str,
    family: Any,
    cohorts: Sequence[Any] = (),
    on_saved: Optional[Callable[[], None]] = None,
) -> FamilyStatusChip:
    """The status chip; call ``.refresh()`` after a curation change."""
    return FamilyStatusChip(data_dir, family_id, family, cohorts, on_saved)


def render_family_status_summary(
    counts: fs.StatusCounts, *, compact: bool = False, mark: str = "status-summary"
) -> None:
    """A stacked bar of the families by status, and the yield under it.

    ``compact`` for the home cards: a thinner bar and a shorter line.
    """
    total = counts.analysable
    with ui.column().classes("w-full min-w-0 gap-1").mark(mark):
        if total:
            with ui.element("div").classes(
                "flex w-full overflow-hidden rounded" + (" h-2" if compact else " h-3")
            ):
                for key, n in counts.by_status.items():
                    if not n:
                        continue
                    info = fs.status_info(key)
                    ui.element("div").style(
                        f"width: {100 * n / total:.2f}%;"
                        f" background-color: {info['color']}"
                    ).tooltip(f"{info['label']}: {n:,}")
        with ui.row().classes("w-full min-w-0 items-center gap-x-3 gap-y-0 wrap"):
            ui.label(counts.yield_text(compact)).classes(
                "text-xs text-gray-700" + (" truncate min-w-0" if compact else "")
            ).mark(f"{mark}-yield").tooltip(
                "Solved and partially solved families, of the families"
                " reviewed (candidate, partially solved, solved, unsolved)."
                " Analysable: every family with an affected member or a"
                " status, but those not analysable."
            )
            if not compact:
                for key, n in counts.by_status.items():
                    if n:
                        info = fs.status_info(key)
                        with ui.row().classes("items-center gap-1 no-wrap"):
                            ui.element("span").classes(
                                "inline-block w-2 h-2 rounded-sm"
                            ).style(f"background-color: {info['color']}")
                            ui.label(f"{info['label']} {n:,}").classes(
                                "text-xs text-gray-600 whitespace-nowrap"
                            )
            if counts.reopened:
                ui.label(f"{counts.reopened:,} reopened").classes(
                    "text-xs text-amber-800 whitespace-nowrap"
                ).tooltip("A pedigree edit after the status: review again")
            if counts.excluded and not compact:
                ui.label(f"{counts.excluded:,} not counted").classes(
                    "text-xs text-gray-400 whitespace-nowrap"
                ).tooltip("Not analysable, or no affected member and no status")
