"""Counts drawn on tabs (``utils/review_counts.py``).

:class:`TabBadge` is the family page's: content that waits for a curator,
computed off the event loop and drawn when it is in, so a page never waits on
it; :meth:`TabBadge.update` recomputes it after a save. A count of 0 draws
nothing.

:class:`TabPills` is the cohort page's: small badges after the tab's label,
drawn from whatever list of :class:`Pill` the page computed last.
"""

import asyncio
import logging
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional

from nicegui import ui

from genetics_viz.utils.review_counts import (
    Pill,
    ReviewCount,
    denovo_to_review,
    qc_family_to_review,
    str_to_review,
    upd_to_review,
)
from genetics_viz.utils.sample_qc import severity_info

logger = logging.getLogger(__name__)


class TabBadge:
    """The badge of one ``ui.tab``, marked ``mark`` for the tests."""

    def __init__(
        self, tab: ui.tab, mark: str, compute: Callable[[], ReviewCount]
    ) -> None:
        self._tab = tab
        self._mark = mark
        self._compute = compute
        self._count: Optional[ReviewCount] = None
        with tab:

            @ui.refreshable
            def draw() -> None:
                self._render()

            draw()
        self._draw = draw

    def _render(self) -> None:
        count = self._count
        if count is None or not count.shown:
            return
        ui.badge(
            "!" if count.held else str(count.n),
            color=severity_info(count.worst)["color"],
        ).props("floating").tooltip(count.tooltip).mark(self._mark)

    async def update(self) -> None:
        """Recompute the count in a thread, then redraw."""
        try:
            self._count = await asyncio.to_thread(self._compute)
        except Exception:
            logger.warning("Failed to count %s", self._mark, exc_info=True)
            self._count = None
        self._draw.refresh()

    def refresh(self) -> None:
        """:meth:`update` from a sync callback (a save, a review closed).

        The timer hangs from the tab, which outlives the dialog that saved.
        """
        with self._tab:
            ui.timer(0.05, self.update, once=True)


class TabPills:
    """Badges after a ``ui.tab``'s label, each marked ``mark`` for the tests.

    Not floating: a tab may carry two (high and medium), and the cohort tabs
    have room for them in line.
    """

    def __init__(self, tab: ui.tab, mark: str) -> None:
        self._mark = mark
        self.pills: List[Pill] = []
        with tab:

            @ui.refreshable
            def draw() -> None:
                self._render()

            draw()
        self._draw = draw

    def _render(self) -> None:
        for pill in self.pills:
            badge = ui.badge(pill.text, color=pill.color).classes(
                "ml-1 text-[11px] leading-none"
            )
            if pill.icon:
                badge.props("rounded").classes("px-[3px]")
                with badge:
                    ui.icon(pill.icon, size="12px")
            if pill.tooltip:
                badge.tooltip(pill.tooltip)
            badge.mark(self._mark)

    def set(self, pills: List[Pill]) -> None:
        """Draw ``pills``, redrawing only when they changed."""
        if pills != self.pills:
            self.pills = list(pills)
            self._draw.refresh()


def family_review_badges(
    data_dir: Path,
    family_id: str,
    tabs: Mapping[str, ui.tab],
    available: Mapping[str, bool],
    show_review: bool,
) -> Dict[str, TabBadge]:
    """The family page's badges, counted off-thread once the page is drawn.

    ``tabs`` and ``available`` are keyed by tab slug (``qc``, ``wombat``,
    ``str``, ``upd``). Sample QC is for everyone, as the header QC badge;
    the others count what only a curator settles, so only ``show_review``
    (``can_write()``, resolved on the event loop) gets them.
    """
    has_upd = bool(available.get("upd"))
    counts: Dict[str, tuple] = {
        "qc": (
            "qc-tab-count",
            lambda: qc_family_to_review(data_dir, family_id, has_upd),
        ),
    }
    if show_review:
        counts.update(
            wombat=("wombat-tab-count", lambda: denovo_to_review(data_dir, family_id)),
            str=("str-tab-count", lambda: str_to_review(data_dir, family_id)),
            upd=("upd-tab-count", lambda: upd_to_review(data_dir, family_id)),
        )
    badges = {
        key: TabBadge(tabs[key], mark, compute)
        for key, (mark, compute) in counts.items()
        if available.get(key) and key in tabs
    }

    async def update_all() -> None:
        for badge in badges.values():
            await badge.update()

    if badges:
        ui.timer(0.2, update_all, once=True)
    return badges
