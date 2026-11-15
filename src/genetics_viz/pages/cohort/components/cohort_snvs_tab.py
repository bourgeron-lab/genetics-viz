"""Cohort-wide SNVs/Indels tab: de novo candidates across families.

The per-variant shortcut the STR, UPD and ROH tabs already give their kinds:
the ``denovo`` step's candidates of the children the filters leave, then the
genes hit in two or more families. Strict HIGH/MODERATE candidates by
default; inherited variants are the Search page's. The cohort tables and the
curation files are read off-thread; the per-child QC is the De novo QC tab's.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from nicegui import ui

from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils import cohort_snvs, denovo
from genetics_viz.utils.export import export_filename, rows_to_tsv
from genetics_viz.utils.gene_scoring import build_gene_badges, get_gene_scorer
from genetics_viz.utils.validation_status import FILTER_STATUSES

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_snvs", "render_cohort_snvs_tab"]

_NUM = {"sortable": True, "sorting": "numerical"}

_CANDIDATE_EXPORT = [
    "FID", "child", "Variant", "var_type", "genes", "consequence", "impact",
    "strict", "flags", "origin", "child_vaf", "recurrence_families", "Validation",
    "Diagnostic",
]  # fmt: skip
_GENE_EXPORT = [
    "gene", "n_families", "n_children", "n_candidates", "families", "consequences",
    "impact",
]  # fmt: skip


def probe_cohort_snvs(cohort: Any) -> bool:
    """The tab is enabled when the cohort has the denovo step's tables."""
    return denovo.probe_cohort_denovo(cohort)


def _keys(filtered_state: Dict[str, Any]) -> Set[Tuple[str, str]]:
    return {
        (str(ind.get("FID")), str(ind["Sample ID"]))
        for ind in filtered_state.get("individuals", [])
        if ind.get("Sample ID")
    }


def _load(cohort: Any, data_dir: Any, scorer: Any) -> Dict[str, Any]:
    """Blocking: every candidate as a table row, with curation and gene badges."""
    data = denovo.load_cohort_denovo(cohort)
    if data is None:
        return {"label": None, "rows": []}
    rows = data.rows()
    families = {str(r.get("family_id")) for r in rows}
    curation = denovo.load_curation(data_dir, families)
    table = cohort_snvs.candidate_rows(rows, None, curation)
    for r in table:
        r["genes_badges"] = build_gene_badges(r["_genes"], scorer)
    return {"label": data.label, "rows": table}


def _candidate_columns(cohort: Any) -> List[Dict[str, Any]]:
    return [
        {
            "id": "FID",
            "header": "Family",
            "cellType": "link",
            "href": f"/cohort/{cohort.name}/family/{{FID}}?tab=snvs&wcfg=denovo",
            "sortable": True,
        },
        {"id": "child", "header": "Child", "sortable": True},
        {"id": "Variant", "header": "Variant", "sortable": True, "sorting": "genomic"},
        {"id": "var_type", "header": "Type", "sortable": True},
        {"id": "genes", "header": "Genes", "cellType": "badge_list", "sortable": True},
        {"id": "consequence", "header": "Consequence", "sortable": True},
        {
            "id": "impact",
            "header": "Impact",
            "cellType": "badge_list",
            "sortable": True,
        },
        {
            "id": "strict",
            "header": "Strict",
            "cellType": "badge_list",
            "sortable": True,
        },
        {"id": "flags", "header": "Flags", "cellType": "badge_list", "sortable": True},
        {"id": "origin", "header": "Origin", "sortable": True},
        {"id": "child_vaf", "header": "VAF", **_NUM},
        {"id": "recurrence_families", "header": "Other families", **_NUM},
        {
            "id": "Validation",
            "header": "Validation",
            "cellType": "badge_list",
            "sortable": True,
        },
        {"id": "Diagnostic", "header": "Diagnostic", "sortable": True},
    ]


def _gene_columns() -> List[Dict[str, Any]]:
    return [
        {"id": "gene", "header": "Gene", "cellType": "badge_list", "sortable": True},
        {"id": "n_families", "header": "Families", **_NUM},
        {"id": "n_children", "header": "Children", **_NUM},
        {"id": "n_candidates", "header": "Candidates", **_NUM},
        {"id": "families", "header": "Family IDs", "sortable": True},
        {"id": "consequences", "header": "Consequences", "sortable": True},
        {"id": "impact", "header": "Worst impact", "cellType": "badge_list"},
    ]


def render_cohort_snvs_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort-wide SNVs/Indels tab."""
    if not probe_cohort_snvs(cohort):
        ui.label(
            "No cohort de novo tables: the SNVs/Indels view lists the denovo"
            " step's candidates"
        ).classes("text-gray-500 italic")
        return
    # The scorer reads the user's storage: take it on the event loop.
    scorer = get_gene_scorer()
    state: Dict[str, Any] = {"loaded": None}
    filters: Dict[str, Any] = {
        "strict": True,
        "impacts": list(cohort_snvs.DEFAULT_IMPACTS),
        "statuses": [],
    }

    def shown() -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """The candidates the filters leave, and their recurrent genes."""
        loaded = state["loaded"] or {"rows": []}
        keys = _keys(filtered_state)
        rows = [r for r in loaded["rows"] if (r["FID"], r["child"]) in keys]
        rows = cohort_snvs.filter_candidates(
            rows,
            strict_only=filters["strict"],
            impacts=filters["impacts"],
            statuses=filters["statuses"],
        )
        genes = cohort_snvs.recurrent_genes(rows)
        for g in genes:
            g["gene_badges"] = build_gene_badges([g["gene"]], scorer)
        return rows, genes

    async def load() -> None:
        try:
            state["loaded"] = await asyncio.to_thread(
                _load, cohort, store.data_dir, scorer
            )
        except Exception:
            logger.warning(
                "Failed to read the de novo candidates of %s",
                cohort.name,
                exc_info=True,
            )
            state["loaded"] = {"label": None, "rows": []}
        render_body.refresh()

    with ui.card().classes("w-full min-w-0"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("De novo SNVs and indels").classes(
                "text-lg font-semibold text-blue-700"
            )
            ui.label(
                "candidates of the children shown (denovo step); the per-child QC"
                " is on the De novo QC tab"
            ).classes("text-xs text-gray-500")
            ui.space()
            ui.button(
                "Search inherited variants",
                icon="search",
                on_click=lambda: ui.navigate.to("/search"),
            ).props("flat dense no-caps color=blue-8").mark("snvs-search")

            def label() -> Optional[str]:
                return (state["loaded"] or {}).get("label")

            download_button(
                [
                    Download(
                        "candidates",
                        "The candidates shown, with their flags and curation",
                        lambda: export_filename(
                            cohort.name, "snvs", "candidates", tag=label()
                        ),
                        lambda: rows_to_tsv(shown()[0], _CANDIDATE_EXPORT),
                    ),
                    Download(
                        "genes",
                        "Genes hit in two or more families",
                        lambda: export_filename(
                            cohort.name, "snvs", "recurrent_genes", tag=label()
                        ),
                        lambda: rows_to_tsv(shown()[1], _GENE_EXPORT),
                    ),
                ],
                mark="download-snvs",
            )

        with ui.row().classes("items-center gap-4 w-full"):

            def on_strict(e: Any) -> None:
                filters["strict"] = bool(e.value)
                render_body.refresh()

            def on_impacts(e: Any) -> None:
                filters["impacts"] = list(e.value or [])
                render_body.refresh()

            def on_statuses(e: Any) -> None:
                filters["statuses"] = list(e.value or [])
                render_body.refresh()

            ui.switch("Strict only", value=True, on_change=on_strict).tooltip(
                "The step's strict candidates: no flag but a sibling sharing it"
            ).mark("snvs-strict")
            ui.select(
                list(cohort_snvs.IMPACTS),
                value=list(cohort_snvs.DEFAULT_IMPACTS),
                multiple=True,
                label="Impact",
                on_change=on_impacts,
            ).props("outlined dense use-chips").classes("w-72").mark("snvs-impacts")
            ui.select(
                list(FILTER_STATUSES),
                value=[],
                multiple=True,
                label="Validation (all)",
                on_change=on_statuses,
            ).props("outlined dense use-chips").classes("w-72").mark("snvs-statuses")

        @ui.refreshable
        def render_body() -> None:
            if state["loaded"] is None:
                with ui.row().classes("items-center gap-2"):
                    ui.spinner(size="sm")
                    ui.label("Reading the cohort's de novo candidates…").classes(
                        "text-xs text-gray-500"
                    )
                return
            rows, genes = shown()
            n_families = len(cohort_snvs.families_of(rows))
            ui.label(
                f"{len(rows):,} candidate{'s' if len(rows) != 1 else ''} in"
                f" {n_families:,} famil{'ies' if n_families != 1 else 'y'};"
                f" {len(genes):,} gene{'s' if len(genes) != 1 else ''} hit in two"
                " or more families."
            ).classes("text-xs text-gray-500").mark("snvs-note")

            ui.label("Genes hit in two or more families").classes(
                "text-base font-semibold text-blue-700 mt-2"
            )
            if genes:
                DataTable(
                    columns=_gene_columns(),
                    rows=genes,
                    row_key="gene",
                    pagination={"rowsPerPage": 10},
                )
            else:
                ui.label("No gene is hit in two families with these filters").classes(
                    "text-sm text-gray-500 italic"
                )

            ui.label("Candidates").classes("text-base font-semibold text-blue-700 mt-4")
            if rows:
                DataTable(
                    columns=_candidate_columns(cohort),
                    rows=rows,
                    row_key="_key",
                    pagination={"rowsPerPage": 20},
                )
            else:
                ui.label("No candidate among the children shown").classes(
                    "text-sm text-gray-500 italic"
                )

        render_body()
        refreshers.append(render_body.refresh)
        ui.timer(0.05, load, once=True)
