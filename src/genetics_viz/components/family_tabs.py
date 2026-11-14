"""The family pages' tab bar: the cohort page's groups, in the same order.

Both family pages draw it (``pages/cohort/family.py`` and
``standalone_family.py``). The order is the curator's: check the family is what
the pedigree says, then its variants, then its genetic background.
"""

from typing import Dict, List, Mapping, Optional, Tuple

from nicegui import ui

from genetics_viz.components.tab_groups import GROUPED_TABS_PROPS, tab_group

#: The family tab rows' props (the analysis tabs and the General tab).
TABS_PROPS = GROUPED_TABS_PROPS

#: The analysis tabs by group, in order: link value (``tab=``), label, tooltip.
FAMILY_TAB_GROUPS: List[Tuple[str, List[Tuple[str, str, str]]]] = [
    (
        "Checks",
        [
            (
                "qc",
                "Sample QC",
                "Relatedness, sex and contamination (somalier); DNA source and library",
            ),
        ],
    ),
    (
        "Variants",
        [
            (
                "snvs",
                "SNVs/Indels",
                "Small variants: de novo candidates, then the Wombat tables",
            ),
            ("svs", "SVs", "Structural variants: WisecondorX and Vedette calls"),
            ("str", "STR", "Short tandem repeats"),
            ("upd", "UPD", "Uniparental disomy"),
        ],
    ),
    (
        "Background",
        [
            ("ancestry", "Ancestry", ""),
            ("roh", "ROH", "Runs of homozygosity"),
            ("pgs", "Polygenic Scores", ""),
        ],
    ),
]

#: The analysis tabs, in order: link value, label, tooltip.
FAMILY_TABS: List[Tuple[str, str, str]] = [
    tab for _caption, tabs in FAMILY_TAB_GROUPS for tab in tabs
]

#: A tab's label by link value; the ``update:model-value`` event reports labels.
TAB_LABELS: Dict[str, str] = {slug: label for slug, label, _tip in FAMILY_TABS}
#: A tab's link value by label.
TAB_SLUGS: Dict[str, str] = {label: slug for slug, label in TAB_LABELS.items()}

#: Link values of tabs since renamed: ``tab=wombat`` opens SNVs/Indels.
TAB_ALIASES: Dict[str, str] = {"wombat": "snvs"}


def canonical_tab(value: Optional[str]) -> str:
    """A linked tab's current link value (an older name is mapped)."""
    value = value or ""
    return TAB_ALIASES.get(value, value)


def family_tab_bar(available: Mapping[str, bool]) -> Tuple[ui.tabs, Dict[str, ui.tab]]:
    """The analysis tab bar, marked ``family-tabs``, and its tabs by link value.

    A tab that ``available`` maps to false is disabled; one it leaves out is not.
    """
    with ui.tabs().props(TABS_PROPS).classes("w-full").mark("family-tabs") as bar:
        tabs: Dict[str, ui.tab] = {}
        for index, (caption, group) in enumerate(FAMILY_TAB_GROUPS):
            with tab_group(caption, first=index == 0):
                for slug, label, tip in group:
                    tab = ui.tab(label)
                    if tip:
                        tab.tooltip(tip)
                    if not available.get(slug, True):
                        tab.props("disable")
                    tabs[slug] = tab
    return bar, tabs
