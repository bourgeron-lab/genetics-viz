"""Cohort-wide STR sub-tab for the cohort page's right-hand panel.

ExpansionHunter only, for now. Follows the filters of the individuals table:

1. per locus, how many of the shown individuals are pathogenic, carrier or
   intermediate, beside the founders' 99.9th percentile and maximum - those
   come from ``loci.tsv`` and describe the whole cohort, not the filter;
2. the shown individuals' flagged genotypes, linked to their family.

A locus row opens the founders' distribution with the flagged individuals
marked on it. The cohort genotype table is several megabytes, so it is read
off-thread before anything is drawn.
"""

import asyncio
import logging
from typing import Any, Callable, Dict, List, Optional, Tuple

from nicegui import ui

from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.str_eh import (
    add_locus_badges,
    CohortContext,
    cohort_flagged_rows,
    cohort_locus_rows,
    find_eh_labels,
    get_cohort_eh_dir,
    load_cohort_context,
    locus_distribution,
    probe_cohort_eh_data,
)
from genetics_viz.utils.str_plots import distribution_option

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_eh_data", "render_cohort_str_tab"]

#: Members marked on a locus plot at most: one row each, so a few hundred
#: flagged individuals (ARX_2) would make the chart metres tall.
_MAX_MARKED = 25


def _table_box() -> ui.element:
    """A block box the table can scroll inside (see cohort_roh_tab)."""
    return ui.element("div").classes("w-full min-w-0")


def _keys(filtered_state: Dict[str, Any]) -> List[Tuple[Optional[str], str]]:
    return [
        (ind.get("FID"), ind["Sample ID"])
        for ind in filtered_state.get("individuals", [])
        if ind.get("Sample ID")
    ]


def render_cohort_str_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort-wide STR sub-tab.

    Args:
        store: DataStore instance
        cohort: Cohort instance
        filtered_state: Dict with 'individuals' key holding the filtered rows
        refreshers: List to append the refresh callback to
    """
    eh_dir = get_cohort_eh_dir(cohort)
    labels = find_eh_labels(eh_dir, cohort.name)
    if not labels:
        ui.label(f"No cohort-level ExpansionHunter results in: {eh_dir}").classes(
            "text-gray-500 italic"
        )
        return

    state: Dict[str, Any] = {"label": labels[0], "context": None, "loading": True}

    async def load() -> None:
        state["loading"] = True
        label = state["label"]
        try:
            context = await asyncio.to_thread(load_cohort_context, cohort, label)
        except Exception:
            logger.warning(
                "Failed to read cohort EH for %s", cohort.name, exc_info=True
            )
            context = None
        if label != state["label"]:
            return  # the catalog changed while reading
        state["context"] = context
        state["loading"] = False
        render_body.refresh()

    with ui.card().classes("w-full min-w-0"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("ExpansionHunter").classes("text-lg font-semibold text-blue-700")
            if len(labels) == 1:
                ui.label(f"catalog {labels[0]}").classes("text-xs text-gray-500")

        if len(labels) > 1:

            def on_label(e: Any) -> None:
                state.update(label=e.value, context=None, loading=True)
                render_body.refresh()
                ui.timer(0.05, load, once=True)

            ui.select(
                options=labels, value=labels[0], label="Catalog", on_change=on_label
            ).props("outlined dense").classes("w-full mb-2")

        @ui.refreshable
        def render_body() -> None:
            context: Optional[CohortContext] = state["context"]
            if state["loading"]:
                with ui.row().classes("items-center gap-2"):
                    ui.spinner(size="sm")
                    ui.label("Reading the cohort's ExpansionHunter tables…").classes(
                        "text-xs text-gray-500"
                    )
                return
            if context is None or context.genotypes is None:
                ui.label(
                    f"No cohort genotype table for catalog {state['label']}"
                ).classes("text-gray-500 italic")
                return

            keys = _keys(filtered_state)
            locus_rows = cohort_locus_rows(context, keys)
            flagged = cohort_flagged_rows(context, keys)
            scorer = get_gene_scorer()
            add_locus_badges(locus_rows, scorer)
            add_locus_badges(flagged, scorer)
            shown = len({k for k in keys})
            ui.label(
                f"{shown} individuals shown, {len(flagged)} flagged genotypes."
                + (
                    f" Founders' statistics: all {context.n_founders} founders."
                    if context.n_founders
                    else ""
                )
            ).classes("text-xs text-gray-500")

            num = {"sortable": True, "sorting": "numerical"}

            def on_locus(e: Dict[str, Any]) -> None:
                locus = (e.get("row") or {}).get("locus")
                if locus:
                    _show_locus_dialog(context, locus, flagged)

            ui.label("Loci").classes("text-base font-semibold text-blue-700 mt-2")
            with _table_box():
                DataTable(
                    columns=[
                        {
                            "id": "actions",
                            "header": "",
                            "cellType": "action",
                            "actionName": "view_locus",
                            "actionIcon": "insights",
                            "actionColor": "#1976d2",
                            "actionTooltip": "Founders' distribution",
                            "sortable": False,
                        },
                        {
                            "id": "locus",
                            "header": "Locus",
                            "sortable": True,
                            "filter": {"type": "text"},
                            "cellType": "gene_badge",
                            "badgesField": "LocusBadges",
                        },
                        {"id": "n_pathogenic", "header": "Path.", **num},
                        {"id": "n_carrier", "header": "Carr.", **num},
                        {"id": "n_intermediate", "header": "Int.", **num},
                        {"id": "founders_p999", "header": "p99.9", **num},
                        {"id": "founders_max", "header": "Max", **num},
                        {
                            "id": "flags",
                            "header": "Flags",
                            "cellType": "badge_list",
                            "sortable": False,
                        },
                    ],  # fmt: skip
                    rows=locus_rows,
                    row_key="id",
                    pagination={"rowsPerPage": 10},
                    on_row_action=on_locus,
                )
            ui.label(
                "Counts: the individuals shown. p99.9 and max: the founders' longest"
                " allele, whole cohort. common: 5% or more of the cohort flagged."
            ).classes("text-xs text-gray-500")

            ui.label("Flagged genotypes").classes(
                "text-base font-semibold text-blue-700 mt-4"
            )
            if not flagged:
                ui.label("No flagged genotype among the individuals shown").classes(
                    "text-gray-500 italic"
                )
                return
            with _table_box():
                DataTable(
                    columns=[
                        {
                            "id": "sample",
                            "header": "Sample",
                            "sortable": True,
                            "filter": {"type": "text"},
                        },
                        {
                            "id": "FID",
                            "header": "Family",
                            "cellType": "link",
                            "href": f"/cohort/{cohort.name}/family/{{FID}}",
                            "sortable": True,
                        },
                        {
                            "id": "locus",
                            "header": "Locus",
                            "sortable": True,
                            "filter": {"type": "text"},
                            "cellType": "gene_badge",
                            "badgesField": "LocusBadges",
                        },
                        {
                            "id": "status",
                            "header": "Status",
                            "cellType": "badge_list",
                            "sortField": "status_rank",
                            **num,
                        },
                        {
                            "id": "alleles",
                            "header": "Alleles",
                            "cellType": "badge_list",
                            "sortField": "longest",
                            **num,
                        },
                        {"id": "pct", "header": "Pct", **num},
                    ],  # fmt: skip
                    rows=flagged,
                    row_key="id",
                    pagination={"rowsPerPage": 10},
                    initial_sorting=[{"id": "status", "desc": True}],
                )

        render_body()
        refreshers.append(render_body.refresh)
        ui.timer(0.05, load, once=True)


def _show_locus_dialog(
    context: CohortContext, locus: str, flagged: List[Dict[str, Any]]
) -> None:
    """The founders' distribution at one locus, the shown flagged individuals on it."""
    info = next(
        (r for r in (context.loci.to_dicts() if context.loci is not None else [])
         if r["locus"] == locus),
        {},
    )  # fmt: skip
    at_locus = sorted(
        (r for r in flagged if r["locus"] == locus),
        key=lambda r: -(r["longest"] or 0),
    )
    marked = [
        (
            f"{r['sample']} ({r['FID']})",
            r["allele_values"],
            r["status"],
        )
        for r in at_locus[:_MAX_MARKED]
    ]  # fmt: skip
    title = f"{locus} in cohort {context.name}"
    option = distribution_option(
        locus_distribution(context.genotypes, locus), info, marked, title=title
    )
    if len(at_locus) > _MAX_MARKED:
        option["title"]["subtext"] += (
            f" Marked: the {_MAX_MARKED} longest of {len(at_locus)} flagged."
        )

    with ui.dialog().props("full-width") as dialog, ui.card().classes("w-full"):
        with ui.row().classes("w-full items-center"):
            ui.label(locus).classes("text-lg font-semibold text-blue-700")
            ui.label(
                f"{info.get('disease') or 'no disease'} · {info.get('inheritance') or '–'}"
            ).classes("text-sm text-gray-600")
            ui.space()
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")
        if not marked:
            ui.label("No flagged individual shown at this locus.").classes(
                "text-sm text-gray-500"
            )
        height = 420 + 22 * len(marked)
        # shrink-0: the dialog card is a flex column and would squeeze the chart.
        ui.echart(option).classes("w-full shrink-0").style(f"height: {height}px")
    dialog.open()
