"""Asking for a ghfc-ngs run: the dialog, the request button, the Workflow list.

Curators and administrators only, and only for a data directory with run
requests turned on (Admin → Directories). The request is a file in the cohort
directory (``utils/run_requests.py``); the ghfc-ngs-swiftbar dashboard shows it
to the person who runs the pipeline, who starts or dismisses it there.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any, Callable, List, Optional, Sequence

from nicegui import ui

from genetics_viz.utils import run_requests as rr
from genetics_viz.utils.auth import can_write, get_current_user, is_admin
from genetics_viz.utils.cohort_state import step_label

logger = logging.getLogger(__name__)

#: Chip colour and label of each request state, and of a started run's outcome.
_STATE = {
    rr.OPEN: ("orange", "waiting"),
    rr.WITHDRAWN: ("grey", "withdrawn"),
    rr.DISMISSED: ("grey-7", "dismissed"),
    rr.STARTED: ("blue", "started"),
}
_OUTCOME = {
    "running": ("blue", "running"),
    "success": ("green", "succeeded"),
    "failed": ("red", "failed"),
    "interrupted": ("red", "interrupted"),
    "died": ("red", "died"),
}


def requests_shown(store: Any) -> bool:
    """Whether this user sees run requests for ``store``'s data directory."""
    return can_write() and rr.requests_enabled(store.data_dir)


def request_pills(waiting: int) -> List[Any]:
    """The Workflow button's badge for requests not yet handled."""
    from genetics_viz.utils.review_counts import Pill

    if not waiting:
        return []
    word = "request" if waiting == 1 else "requests"
    return [Pill(str(waiting), "orange", f"{waiting} run {word} waiting")]


def _when(moment: Any) -> str:
    return moment.strftime("%Y-%m-%d %H:%M") if moment is not None else "?"


def _steps_words(steps: Optional[Sequence[str]]) -> str:
    if steps is None:
        return "normal run"
    return ", ".join(step_label(s) for s in steps)


def open_request_dialog(
    cohort: Any,
    *,
    steps: Optional[Sequence[str]] = None,
    clean: bool = False,
    families: Sequence[str] = (),
    origin: str = "",
    reason: str = "",
    on_saved: Optional[Callable[[], Any]] = None,
) -> None:
    """The "Request a run" dialog, prefilled. ``steps=None`` starts on a
    normal run. Saving checks the role again, in the writer."""
    allowed = rr.allowed_steps(cohort)
    chosen = {s: (steps is not None and s in steps) for s in allowed}
    form = {"normal": steps is None or not any(chosen.values()), "clean": clean}

    with ui.dialog() as dialog, ui.card().classes("w-[34rem] max-w-full gap-3"):
        dialog.mark("run-request-dialog")
        ui.label(f"Request a ghfc-ngs run of {cohort.name}").classes(
            "text-lg font-semibold"
        )
        ui.label(
            "The person who runs the pipeline sees the request in the ghfc-ngs"
            " dashboard, and starts it (possibly with other steps) or dismisses it."
            " Each step only redoes what is out of date."
        ).classes("text-xs text-gray-600")
        if families:
            ui.label("Families: " + ", ".join(families)).classes(
                "text-xs text-gray-600"
            )

        mode = (
            ui.radio(
                {
                    "steps": "These steps",
                    "normal": "A normal run (every step of the params file)",
                },  # fmt: skip
                value="normal" if form["normal"] else "steps",
            )
            .props("dense")
            .classes("text-sm")
            .mark("run-request-mode")
        )
        boxes = ui.row().classes("w-full gap-x-4 gap-y-0 flex-wrap pl-6")
        with boxes:
            for step in allowed:

                def make_toggle(name: str):
                    def toggle(e: Any) -> None:
                        chosen[name] = bool(e.value)
                        update()

                    return toggle

                ui.checkbox(
                    step_label(step), value=chosen[step], on_change=make_toggle(step)
                ).props("dense").classes("text-sm").mark(f"run-request-step-{step}")

        def set_clean(e: Any) -> None:
            form["clean"] = bool(e.value)
            update()

        ui.checkbox("--clean-stale-families", value=clean, on_change=set_clean).props(
            "dense"
        ).classes("text-sm font-mono").mark("run-request-clean").tooltip(
            "For families whose sequenced members changed: deletes their joint"
            " call and what follows it, so they are called again in full."
        )
        reason_box = (
            ui.textarea("Why (shown with the request)", value=reason)
            .props("outlined dense autogrow")
            .classes("w-full")
            .mark("run-request-reason")
        )
        command = ui.label("").classes(
            "text-xs font-mono bg-gray-100 rounded px-2 py-1 w-full break-all"
        )
        command.mark("run-request-command")

        def selected() -> Optional[List[str]]:
            if mode.value == "normal":
                return None
            return [s for s in allowed if chosen[s]]

        def update() -> None:
            boxes.set_visibility(mode.value == "steps")
            picked = selected()
            if picked is not None and not picked:
                command.set_text("Choose at least one step")
            else:
                command.set_text(rr.command_line(cohort.name, picked, form["clean"]))

        mode.on_value_change(lambda _: update())
        update()

        async def send() -> None:
            try:
                request = await rr.write_request_async(
                    cohort,
                    selected(),
                    clean=form["clean"],
                    reason=reason_box.value or "",
                    user=get_current_user(),
                    families=list(families),
                    origin=origin,
                )
            except (ValueError, PermissionError) as e:
                ui.notify(str(e), type="negative")
                return
            except Exception:
                logger.exception("Failed to write a run request for %s", cohort.name)
                ui.notify("The request could not be written", type="negative")
                return
            ui.notify(f"Run requested: {_steps_words(request.steps)}", type="positive")
            dialog.close()
            if on_saved is not None:
                result = on_saved()
                if asyncio.iscoroutine(result):
                    await result

        with ui.row().classes("w-full justify-end gap-2"):
            ui.button("Cancel", on_click=dialog.close).props("flat")
            ui.button("Send request", icon="send", on_click=send).props(
                "unelevated color=primary"
            ).mark("run-request-send")
    dialog.open()


async def _withdraw(cohort: Any, request_id: str, after: Callable[[], Any]) -> None:
    try:
        await rr.withdraw_request_async(
            cohort.pedigree_file.parent, request_id, user=get_current_user()
        )
    except (ValueError, PermissionError) as e:
        ui.notify(str(e), type="negative")
        return
    ui.notify("Request withdrawn", type="info")
    result = after()
    if asyncio.iscoroutine(result):
        await result


def render_request_control(
    cohort: Any,
    *,
    label: str = "Request a run",
    mark: str = "run-request",
    steps: Optional[Sequence[str]] = None,
    clean: bool = False,
    families: Sequence[str] = (),
    origin: str = "",
    reason: str = "",
    on_changed: Optional[Callable[[], Any]] = None,
) -> Callable[[], Any]:
    """A button opening the dialog, or, while the user has an open request for
    the cohort, that request with Withdraw. Returns its reload function."""
    box = ui.row().classes("items-center gap-2 flex-wrap")
    user = get_current_user()

    async def reload() -> None:
        try:
            requests = await asyncio.to_thread(
                rr.read_requests, cohort.pedigree_file.parent
            )
        except Exception:
            logger.warning(
                "Failed to read the run requests of %s", cohort.name, exc_info=True
            )
            requests = []
        mine = rr.open_request_of(requests, user)
        box.clear()
        with box:
            if mine is not None:
                ui.icon("schedule_send").classes("text-orange-700")
                ui.label(
                    f"Run requested {_when(mine.requested_at)}: {_steps_words(mine.steps)}"
                ).classes("text-xs text-orange-900").mark(f"{mark}-open")
                ui.button(
                    "Withdraw", on_click=lambda: _withdraw(cohort, mine.id, changed)
                ).props("flat dense size=sm color=grey-8").mark(
                    f"run-request-withdraw-{mine.id}"
                )
            else:
                ui.button(
                    label,
                    icon="send",
                    on_click=lambda: open_request_dialog(
                        cohort,
                        steps=steps,
                        clean=clean,
                        families=families,
                        origin=origin,
                        reason=reason,
                        on_saved=changed,
                    ),
                ).props("outline dense size=sm color=primary").mark(mark)

    async def changed() -> None:
        await reload()
        if on_changed is not None:
            result = on_changed()
            if asyncio.iscoroutine(result):
                await result

    ui.timer(0.05, reload, once=True)
    return reload


def render_workflow_requests(cohort: Any) -> None:
    """The Workflow tab's run requests: the button, then every curator's
    requests of the last 30 days with what became of them."""
    with (
        ui.column()
        .classes("w-full min-w-0 gap-2 border border-gray-200 rounded-md p-3")
        .mark("workflow-requests")
    ):
        header = ui.row().classes("w-full items-center justify-between gap-2")
        with header:
            ui.label("Run requests").classes("text-sm font-semibold text-gray-700")
        listing = ui.column().classes("w-full min-w-0 gap-1")

        async def load_list() -> None:
            try:
                requests = await asyncio.to_thread(
                    rr.read_requests, cohort.pedigree_file.parent
                )
            except Exception:
                logger.warning("Failed to read the run requests", exc_info=True)
                return
            shown = rr.recent(requests)
            admin, user = is_admin(), get_current_user()
            listing.clear()
            with listing:
                if not shown:
                    ui.label("No request in the last 30 days.").classes(
                        "text-xs text-gray-500 italic"
                    )
                for r in shown:
                    _request_row(cohort, r, admin or r.requested_by == user, refresh)

        async def refresh() -> None:
            await control_reload()
            await load_list()

        with header:
            control_reload = render_request_control(
                cohort, origin="Workflow tab", on_changed=load_list
            )
        ui.timer(0.1, load_list, once=True)


def _request_row(
    cohort: Any, r: rr.RunRequest, may_withdraw: bool, after: Callable[[], Any]
) -> None:
    color, word = _STATE[r.state]
    with (
        ui.column()
        .classes("w-full min-w-0 gap-0 border-t border-gray-100 pt-1")
        .mark(f"run-request-{r.id}")
    ):
        with ui.row().classes("w-full items-center gap-2 no-wrap min-w-0"):
            ui.badge(word, color=color).classes("text-[10px]")
            if r.outcome in _OUTCOME:
                oc, ow = _OUTCOME[r.outcome]
                ui.badge(ow, color=oc).classes("text-[10px]")
            ui.label(
                _steps_words(r.steps) + (" + clean" if r.clean_stale_families else "")
            ).classes("text-xs text-gray-800 truncate min-w-0")
            ui.label(f"{r.requested_by}, {_when(r.requested_at)}").classes(
                "text-xs text-gray-500 whitespace-nowrap"
            )
            if r.is_open and may_withdraw:
                ui.button(
                    "Withdraw", on_click=lambda: _withdraw(cohort, r.id, after)
                ).props("flat dense size=sm color=grey-8").mark(
                    f"run-request-list-withdraw-{r.id}"
                )
        details = [r.reason]
        if r.families:
            details.append("families " + ", ".join(r.families))
        if r.origin:
            details.append(r.origin)
        ui.label(" · ".join(d for d in details if d)).classes(
            "text-xs text-gray-600 break-words"
        )
        h = r.handled
        if r.state == rr.WITHDRAWN:
            ui.label(f"Withdrawn by {r.withdrawn_by}, {_when(r.withdrawn_at)}").classes(
                "text-xs text-gray-500"
            )
        elif h is not None and h.status == "dismissed":
            ui.label(
                f"Dismissed by {h.user}, {_when(h.at)}"
                + (f": {h.note}" if h.note else "")
            ).classes("text-xs text-gray-600")
        elif h is not None:
            text = f"Started by {h.user}, {_when(h.at)}"
            if h.job_id:
                text += f", Slurm job {h.job_id}"
            if h.covered_by:
                text += f", with request {h.covered_by}"
            if r.edited:
                text += f"; sent: {_steps_words(h.steps_sent)}" + (
                    " + clean" if h.clean_sent else ""
                )
            ui.label(text).classes("text-xs text-gray-600")
