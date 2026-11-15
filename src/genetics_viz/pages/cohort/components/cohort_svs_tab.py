"""Cohort-wide SVs tab: large copy-number calls across families.

The cohort's WisecondorX calls of the individuals shown, robust and 1 Mb or
more by default; whole-chromosome gains and losses (the family report's
rule); and regions where calls of two or more families overlap. The table is
scanned lazily and off-thread (``utils/cohort_svs.py``). Vedette has no
cohort table yet.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any, Callable, Dict, List, Set, Tuple

from nicegui import ui

from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import load_validation_map
from genetics_viz.utils import cohort_svs
from genetics_viz.utils.export import export_filename, rows_to_tsv
from genetics_viz.utils.family_report import CONFIG as REPORT_CONFIG
from genetics_viz.utils.gene_scoring import build_gene_badges, get_gene_scorer
from genetics_viz.utils.sample_counts import cohort_wcx_file

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_svs", "render_cohort_svs_tab"]

_NUM = {"sortable": True, "sorting": "numerical"}

_CALL_EXPORT = [
    "FID", "sample", "locus", "size_mb", "call", "ratio", "zscore", "exonic_symbol",
    "Validation",
]  # fmt: skip
_ANEUPLOIDY_EXPORT = ["FID", "sample", "events"]
_REGION_EXPORT = [
    "region", "direction", "size_mb", "n_families", "n_samples", "families",
    "exonic_symbol",
]  # fmt: skip


def probe_cohort_svs(store: Any, cohort: Any) -> bool:
    return cohort_wcx_file(store.data_dir, cohort.name) is not None


def _keys(filtered_state: Dict[str, Any]) -> Set[Tuple[str, str]]:
    return {
        (str(ind.get("FID")), str(ind["Sample ID"]))
        for ind in filtered_state.get("individuals", [])
        if ind.get("Sample ID")
    }


def _load(store: Any, cohort: Any, scorer: Any) -> Dict[str, Any]:
    """Blocking: the calls past the permissive thresholds, as table rows."""
    path = cohort_wcx_file(store.data_dir, cohort.name)
    if path is None:
        return {"rows": [], "file": None}
    calls = cohort_svs.scan_calls(path)
    families_of: Dict[str, List[str]] = {}
    for fid, family in cohort.families.items():
        for sample in family.samples:
            families_of.setdefault(sample.sample_id, []).append(fid)
    validation_file = store.data_dir / "validations" / "svs.tsv"
    rows = cohort_svs.call_rows(
        calls,
        families_of,
        load_validation_map(validation_file) if validation_file.exists() else {},
    )
    for r in rows:
        r["exonic_symbol_badges"] = build_gene_badges(r["_genes"], scorer)
    return {"rows": rows, "file": path.name}


def _call_columns(cohort: Any) -> List[Dict[str, Any]]:
    return [
        {
            "id": "FID",
            "header": "Family",
            "cellType": "link",
            "href": f"/cohort/{cohort.name}/family/{{FID}}?tab=svs",
            "sortable": True,
        },
        {"id": "sample", "header": "Sample", "sortable": True},
        {"id": "locus", "header": "Locus", "sortable": True, "sorting": "genomic"},
        {"id": "size_mb", "header": "Size (Mb)", **_NUM},
        {"id": "call", "header": "Call", "cellType": "badge_list", "sortable": True},
        {"id": "ratio", "header": "log2 ratio", **_NUM},
        {"id": "zscore", "header": "z", **_NUM},
        {"id": "exonic_symbol", "header": "Exonic genes", "cellType": "badge_list"},
        {
            "id": "Validation",
            "header": "Validation",
            "cellType": "badge_list",
            "sortable": True,
        },
    ]


def _region_columns() -> List[Dict[str, Any]]:
    return [
        {"id": "region", "header": "Region", "sortable": True, "sorting": "genomic"},
        {"id": "direction", "header": "Gain / loss", "sortable": True},
        {"id": "size_mb", "header": "Size (Mb)", **_NUM},
        {"id": "n_families", "header": "Families", **_NUM},
        {"id": "n_samples", "header": "Individuals", **_NUM},
        {"id": "families", "header": "Family IDs", "sortable": True},
        {"id": "exonic_symbol", "header": "Exonic genes", "cellType": "badge_list"},
    ]


def render_cohort_svs_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort-wide SVs tab."""
    if not probe_cohort_svs(store, cohort):
        ui.label("No cohort WisecondorX table for this cohort").classes(
            "text-gray-500 italic"
        )
        return
    scorer = get_gene_scorer()
    min_fraction = float(REPORT_CONFIG.get("aneuploidy_min_fraction", 0.75))
    state: Dict[str, Any] = {"loaded": None}
    filters: Dict[str, Any] = {
        "min_kb": cohort_svs.DEFAULT_MIN_KB,
        "permissive": False,
    }

    def shown() -> Tuple[List[Dict], List[Dict], List[Dict]]:
        """The calls the filters leave, their aneuploidies and recurrent regions."""
        rows = (state["loaded"] or {"rows": []})["rows"]
        calls = cohort_svs.filter_calls(
            rows,
            _keys(filtered_state),
            min_kb=filters["min_kb"],
            strengths=cohort_svs.STRENGTHS if filters["permissive"] else ("robust",),
        )
        regions = cohort_svs.recurrent_regions(calls)
        for g in regions:
            g["exonic_symbol_badges"] = build_gene_badges(g["_genes"], scorer)
        return calls, cohort_svs.aneuploidy_rows(calls, min_fraction), regions

    async def load() -> None:
        try:
            state["loaded"] = await asyncio.to_thread(_load, store, cohort, scorer)
        except Exception as e:
            logger.warning("Failed to read %s's SV calls", cohort.name, exc_info=True)
            state["loaded"] = {"rows": [], "file": None, "error": str(e)}
        render_body.refresh()

    with ui.card().classes("w-full min-w-0"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("Copy-number calls").classes("text-lg font-semibold text-blue-700")
            ui.label("WisecondorX, the individuals shown").classes(
                "text-xs text-gray-500"
            )
            ui.space()
            download_button(
                [
                    Download(
                        "calls",
                        "The calls shown",
                        lambda: export_filename(cohort.name, "svs", "calls"),
                        lambda: rows_to_tsv(shown()[0], _CALL_EXPORT),
                    ),
                    Download(
                        "aneuploidies",
                        "Whole-chromosome gains and losses",
                        lambda: export_filename(cohort.name, "svs", "aneuploidies"),
                        lambda: rows_to_tsv(shown()[1], _ANEUPLOIDY_EXPORT),
                    ),
                    Download(
                        "regions",
                        "Regions hit in two or more families",
                        lambda: export_filename(
                            cohort.name, "svs", "recurrent_regions"
                        ),
                        lambda: rows_to_tsv(shown()[2], _REGION_EXPORT),
                    ),
                ],
                mark="download-svs",
            )

        with ui.row().classes("items-center gap-4 w-full"):

            def on_size(e: Any) -> None:
                filters["min_kb"] = int(e.value)
                render_body.refresh()

            def on_permissive(e: Any) -> None:
                filters["permissive"] = bool(e.value)
                render_body.refresh()

            ui.select(
                {kb: f"≥ {kb / 1000:g} Mb" for kb in cohort_svs.SIZE_OPTIONS_KB},
                value=cohort_svs.DEFAULT_MIN_KB,
                label="Size",
                on_change=on_size,
            ).props("outlined dense").classes("w-40").mark("svs-size")
            ui.switch(
                "Include permissive calls", value=False, on_change=on_permissive
            ).tooltip("Robust calls only by default, as the family SVs tab").mark(
                "svs-permissive"
            )

        @ui.refreshable
        def render_body() -> None:
            loaded = state["loaded"]
            if loaded is None:
                with ui.row().classes("items-center gap-2"):
                    ui.spinner(size="sm")
                    ui.label("Reading the cohort's copy-number calls…").classes(
                        "text-xs text-gray-500"
                    )
                return
            if loaded.get("error"):
                ui.label(f"Could not read the calls: {loaded['error']}").classes(
                    "text-sm text-red-700"
                )
                return
            calls, aneuploidies, regions = shown()
            individuals = {(r["FID"], r["sample"]) for r in calls}
            ui.label(
                f"{len(calls):,} call{'s' if len(calls) != 1 else ''} in"
                f" {len(individuals):,} individual{'s' if len(individuals) != 1 else ''}"
                f" ({loaded.get('file')})."
            ).classes("text-xs text-gray-500").mark("svs-note")

            ui.label("Whole-chromosome gains and losses").classes(
                "text-base font-semibold text-blue-700 mt-2"
            )
            if aneuploidies:
                DataTable(
                    columns=[
                        {
                            "id": "FID",
                            "header": "Family",
                            "cellType": "link",
                            "href": f"/cohort/{cohort.name}/family/{{FID}}?tab=svs",
                            "sortable": True,
                        },
                        {"id": "sample", "header": "Sample", "sortable": True},
                        {"id": "events", "header": "Gains / losses"},
                    ],
                    rows=aneuploidies,
                    row_key="_key",
                    pagination={"rowsPerPage": 10},
                )
            else:
                ui.label(
                    f"None: no gain or loss covers {min_fraction * 100:.0f}% of a"
                    " chromosome among the calls shown"
                ).classes("text-sm text-gray-500 italic")

            ui.label("Regions hit in two or more families").classes(
                "text-base font-semibold text-blue-700 mt-4"
            )
            if regions:
                DataTable(
                    columns=_region_columns(),
                    rows=regions,
                    row_key="_key",
                    pagination={"rowsPerPage": 10},
                )
            else:
                ui.label("No region is hit in two families with these filters").classes(
                    "text-sm text-gray-500 italic"
                )

            ui.label("Calls").classes("text-base font-semibold text-blue-700 mt-4")
            if calls:
                DataTable(
                    columns=_call_columns(cohort),
                    rows=calls,
                    row_key="_key",
                    pagination={"rowsPerPage": 20},
                )
            else:
                ui.label("No call among the individuals shown").classes(
                    "text-sm text-gray-500 italic"
                )

        render_body()
        refreshers.append(render_body.refresh)
        ui.timer(0.05, load, once=True)
