"""Tabs in captioned groups, in one row: the cohort and family pages' tab bars.

The tabs carry no icon: a group caption ("Checks", "Variants", "Background")
says what a tab is about better than an icon for "ROH" or "STR" could.

A ``ui.tab`` registers with the nearest ``q-tabs`` above it, whatever lies in
between, so a group is a plain column: its caption over a row of tabs.
"""

from contextlib import contextmanager
from typing import Iterator

from nicegui import ui

#: Mixed-case labels, tabs from the left, 36 px tabs (the caption sits above).
#: ``inline-label`` keeps a tab's content in one line: without it Quasar
#: stacks the label and the cohort page's pills.
GROUPED_TABS_PROPS = "inline-label no-caps dense align=left"


@contextmanager
def tab_group(caption: str, first: bool = False) -> Iterator[None]:
    """A caption over the tabs created in the block, ruled off from the group before.

    The caption is indented as Quasar's 16 px tab padding, so it starts above
    the first label.
    """
    rule = "" if first else " border-l border-gray-200"
    with ui.element("div").classes(f"flex flex-col shrink-0{rule}"):
        ui.label(caption).classes(
            "px-4 pt-1 text-[11px] leading-none text-gray-500 select-none"
        )
        with ui.element("div").classes("flex flex-nowrap"):
            yield
