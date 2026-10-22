"""Wombat tab filters: one config's filter state, its panel, and applying it.

The menus are the shared ones of ``components/filters.py`` and the row tests
those of ``utils/variant_rows.py``, as on the Search page.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List

from nicegui import ui

from genetics_viz.components.page_link import (
    PageLink,
    decode_consequences,
    encode_consequences,
)
from genetics_viz.components.filters import (
    create_consequence_filter_menu,
    create_geneset_filter_menu,
    create_validation_filter_menu,
)
from genetics_viz.utils.genesets import cell_matches, selected_symbols
from genetics_viz.utils.validation_status import FILTER_STATUSES, filter_by_status
from genetics_viz.utils.variant_rows import (
    consequences_selected,
    gnomad_filtered,
    in_lcr,
    likely_mnv,
)
from genetics_viz.utils.vep import VEP_CONSEQUENCES


def _flag(value: bool) -> Dict[str, bool]:
    return {"value": value}


@dataclass
class WombatFilters:
    """One config's filters. ``{"value": ...}`` dicts, as the menus take them.

    The exclusions are on by default; every consequence and validation
    status is selected, and no geneset.
    """

    exclude_lcr: Dict[str, bool] = field(default_factory=lambda: _flag(True))
    exclude_gnomad: Dict[str, bool] = field(default_factory=lambda: _flag(True))
    exclude_mnv: Dict[str, bool] = field(default_factory=lambda: _flag(True))
    genesets: Dict[str, List[str]] = field(default_factory=lambda: {"value": []})
    impacts: Dict[str, List[str]] = field(
        default_factory=lambda: {"value": list(VEP_CONSEQUENCES.keys())}
    )
    validations: Dict[str, List[str]] = field(
        default_factory=lambda: {"value": list(FILTER_STATUSES)}
    )


#: Link keys of the shown config's filters, prefixed ``w`` (Wombat tab).
LINK_KEYS = ("wlcr", "wgnomad", "wmnv", "wgs", "wimp", "wval")


def filters_from_link(link: PageLink) -> WombatFilters:
    """The filters a shared link sets; defaults for what it leaves out."""
    filters = WombatFilters()
    filters.exclude_lcr["value"] = link.get_bool("wlcr", True)
    filters.exclude_gnomad["value"] = link.get_bool("wgnomad", True)
    filters.exclude_mnv["value"] = link.get_bool("wmnv", True)
    if (genesets := link.get_list("wgs")) is not None:
        filters.genesets["value"] = genesets
    if (impacts := link.get_list("wimp")) is not None:
        filters.impacts["value"] = decode_consequences(impacts)
    if (statuses := link.get_list("wval")) is not None:
        filters.validations["value"] = [s for s in statuses if s in FILTER_STATUSES]
    return filters


def filters_to_link(filters: WombatFilters, link: PageLink) -> None:
    """Write the filters into the link, leaving defaults out."""
    link.set("wlcr", filters.exclude_lcr["value"], default=True)
    link.set("wgnomad", filters.exclude_gnomad["value"], default=True)
    link.set("wmnv", filters.exclude_mnv["value"], default=True)
    link.set("wgs", filters.genesets["value"], default=[])
    impacts = filters.impacts["value"]
    link.set(
        "wimp",
        None if set(impacts) == set(VEP_CONSEQUENCES) else encode_consequences(impacts),
    )
    link.set("wval", filters.validations["value"], default=FILTER_STATUSES)


def apply_filters(
    rows: List[Dict[str, Any]], filters: WombatFilters, available_genesets: List[Any]
) -> List[Dict[str, Any]]:
    """The rows that pass every filter."""
    if filters.exclude_lcr["value"]:
        rows = [r for r in rows if not in_lcr(r)]
    if filters.exclude_gnomad["value"]:
        rows = [r for r in rows if not gnomad_filtered(r)]
    if filters.exclude_mnv["value"]:
        rows = [r for r in rows if not likely_mnv(r)]
    if filters.genesets["value"]:
        wanted = selected_symbols(available_genesets, filters.genesets["value"])
        rows = [r for r in rows if cell_matches(r.get("VEP_SYMBOL"), wanted)]
    rows = [r for r in rows if consequences_selected(r, filters.impacts["value"])]
    return filter_by_status(rows, filters.validations["value"])


def render_filter_panel(
    filters: WombatFilters,
    available_genesets: List[Any],
    on_change: Callable[[], Any],
) -> None:
    """The collapsible Filters card: exclusion checkboxes, then the menus."""

    def checkbox(label: str, state: Dict[str, bool]) -> None:
        def changed(e) -> None:
            state["value"] = e.value
            on_change()

        ui.checkbox(label, value=state["value"], on_change=changed)

    with ui.card().classes("w-full p-0 mt-2"):
        with (
            ui.expansion("Filters", icon="filter_list", value=False)
            .classes("w-full")
            .props("header-class='text-lg font-semibold text-blue-700'")
        ):
            with ui.row().classes("gap-6 items-start flex-wrap p-2"):
                # Checkbox exclude filters
                with ui.column().classes("gap-2"):
                    checkbox("Exclude LCR", filters.exclude_lcr)
                    checkbox("Exclude gnomAD filtered", filters.exclude_gnomad)
                    checkbox(
                        "Exclude FP due to MNV (mnv_proba > 0.5)", filters.exclude_mnv
                    )

                # Dropdown menu filters
                with ui.row().classes("gap-2 items-center flex-wrap"):
                    create_geneset_filter_menu(
                        filters.genesets, available_genesets, on_change=on_change
                    )
                    create_consequence_filter_menu(filters.impacts, on_change=on_change)
                    create_validation_filter_menu(
                        all_statuses=list(FILTER_STATUSES),
                        filter_state=filters.validations,
                        on_change=on_change,
                        label="Validation",
                        button_classes="",
                    )
