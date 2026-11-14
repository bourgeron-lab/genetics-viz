"""The family pages' tab bar: an icon on every tab, as on the cohort page.

Both family pages draw it (``pages/cohort/family.py`` and
``standalone_family.py``). A tab the cohort page also has takes its icon.
"""

from typing import Dict, List, Mapping, Optional, Tuple

from nicegui import ui

#: Icon left of the label, mixed-case labels, tabs from the left: the cohort
#: page's tab row. Eight tabs take about 870 px of the 1,185 px a 1,280 px
#: window leaves, so Quasar's 16 px tab padding stays (the cohort page cuts it
#: to 12 for nine tabs). The review badges float in that padding.
TABS_PROPS = "inline-label no-caps align=left"

#: The analysis tabs, in order: link value (``tab=``), label, icon, tooltip.
FAMILY_TABS: List[Tuple[str, str, str, str]] = [
    (
        "snvs",
        "SNVs/Indels",
        "grain",
        "Small variants: de novo candidates, then the Wombat tables",
    ),
    ("svs", "SVs", "call_split", ""),
    ("ancestry", "Ancestry", "public", ""),
    ("pgs", "Polygenic Scores", "functions", ""),
    ("roh", "ROH", "linear_scale", "Runs of homozygosity"),
    ("upd", "UPD", "content_copy", "Uniparental disomy"),
    ("str", "STR", "repeat", "Short tandem repeats"),
    (
        "qc",
        "Sample QC",
        "fact_check",
        "Relatedness, sex and contamination checks (somalier)",
    ),
]

#: A tab's label by link value; the ``update:model-value`` event reports labels.
TAB_LABELS: Dict[str, str] = {slug: label for slug, label, _icon, _tip in FAMILY_TABS}
#: A tab's link value by label.
TAB_SLUGS: Dict[str, str] = {label: slug for slug, label in TAB_LABELS.items()}

#: Link values of tabs since renamed: ``tab=wombat`` opens SNVs/Indels.
TAB_ALIASES: Dict[str, str] = {"wombat": "snvs"}


def canonical_tab(value: Optional[str]) -> str:
    """A linked tab's current link value (an older name is mapped)."""
    value = value or ""
    return TAB_ALIASES.get(value, value)


#: The cohort family page's General tab (notes and diagnostics).
GENERAL_ICON = "info"


def family_tab_bar(available: Mapping[str, bool]) -> Tuple[ui.tabs, Dict[str, ui.tab]]:
    """The analysis tab bar, marked ``family-tabs``, and its tabs by link value.

    A tab that ``available`` maps to false is disabled; one it leaves out is not.
    """
    with ui.tabs().props(TABS_PROPS).classes("w-full").mark("family-tabs") as bar:
        tabs: Dict[str, ui.tab] = {}
        for slug, label, icon, tip in FAMILY_TABS:
            tab = ui.tab(label, icon=icon)
            if tip:
                tab.tooltip(tip)
            if not available.get(slug, True):
                tab.props("disable")
            tabs[slug] = tab
    return bar, tabs
