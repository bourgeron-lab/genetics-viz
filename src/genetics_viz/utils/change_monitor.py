"""Background polling for filesystem changes in data directories."""

from __future__ import annotations

import asyncio
import logging
from collections.abc import Awaitable, Callable

from genetics_viz.models import ChangeReport, DataStore
from genetics_viz.utils.data import get_all_data_stores
from genetics_viz.utils.sharding import clear_sharding_cache

logger = logging.getLogger(__name__)

_poll_task: asyncio.Task[None] | None = None
_change_callbacks: list[Callable[[ChangeReport], Awaitable[None]]] = []


def subscribe(
    callback: Callable[[ChangeReport], Awaitable[None]],
) -> Callable[[], None]:
    """Register a callback for change reports. Returns an unsubscribe function."""
    _change_callbacks.append(callback)

    def unsub() -> None:
        try:
            _change_callbacks.remove(callback)
        except ValueError:
            pass

    return unsub


async def _check_store(store: DataStore, path_str: str) -> ChangeReport | None:
    """Check a single store for changes. Returns a report if changes found."""
    old = store._snapshot
    if old is None:
        return None

    new = await asyncio.to_thread(store.take_snapshot)
    report = DataStore.compare_snapshot(old, new)
    report.data_dir = path_str

    if not report.has_changes:
        return None

    logger.info("Changes detected in %s: %s", path_str, report.summary_lines())
    await asyncio.to_thread(store.reload)
    clear_sharding_cache()
    return report


async def check_now() -> list[ChangeReport]:
    """Run an immediate check on all stores. Returns list of change reports."""
    reports: list[ChangeReport] = []
    for path_str, store in get_all_data_stores().items():
        try:
            report = await _check_store(store, path_str)
            if report is not None:
                reports.append(report)
                await _notify_subscribers(report)
        except Exception:
            logger.exception("Error checking store %s", path_str)
    return reports


async def _notify_subscribers(report: ChangeReport) -> None:
    for cb in list(_change_callbacks):
        try:
            await cb(report)
        except Exception:
            logger.exception("Error in change callback")


async def _poll_loop(interval: float) -> None:
    while True:
        await asyncio.sleep(interval)
        for path_str, store in get_all_data_stores().items():
            try:
                report = await _check_store(store, path_str)
                if report is not None:
                    await _notify_subscribers(report)
            except Exception:
                logger.exception("Error in poll loop for %s", path_str)


def start_polling(interval: float = 30.0) -> None:
    """Start the background polling task. Safe to call multiple times."""
    global _poll_task
    if _poll_task is not None and not _poll_task.done():
        return
    _poll_task = asyncio.get_event_loop().create_task(_poll_loop(interval))
    logger.info("Change monitor started (interval=%ss)", interval)


def stop_polling() -> None:
    """Cancel the background polling task."""
    global _poll_task
    if _poll_task is not None:
        _poll_task.cancel()
        _poll_task = None
        logger.info("Change monitor stopped")
