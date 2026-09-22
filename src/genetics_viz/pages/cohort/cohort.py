"""Cohort detail page - displays pedigree members with filters."""

import asyncio
import csv
from collections import defaultdict
from typing import Any, Callable, Dict, List

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.pages.cohort.components.cohort_ancestry_tab import (
    probe_cohort_ancestry_data,
    render_cohort_ancestry_tab,
)
from genetics_viz.pages.cohort.components.cohort_pgs_tab import (
    probe_cohort_pgs_data,
    render_cohort_pgs_tab,
)
from genetics_viz.pages.cohort.components.docs_panel import render_params_doc_panel
from genetics_viz.pages.cohort.components.params_panel import render_params_panel
from genetics_viz.pages.cohort.components.stats_panel import render_stats_panel
from genetics_viz.pages.cohort.components.status_panel import (
    render_status_panel,
    status_headline,
)
from genetics_viz.utils.auth import check_auth
from genetics_viz.utils.cohort_state import read_cohort_state
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.ghfc_ngs_docs import fetch_params_doc

# Diagnostic priority: higher value wins
_DIAG_PRIORITY = {"pathogenic": 3, "uncertain": 2, "benign": 1}


def _loading(message: str) -> None:
    """Placeholder shown by a lazy tab panel until its content is loaded."""
    with ui.column().classes("w-full items-center justify-center py-16"):
        ui.spinner(size="xl", color="blue")
        ui.label(message).classes("text-lg text-gray-600 mt-4")


def _render_workflow_tab(cohort: Any) -> None:
    """Pipeline status, then the cohort's parameters beside their reference.

    The status panel is open on arrival: it is the reason to come to this tab.
    The two panels below it each scroll internally rather than letting the page
    grow — the reference alone is 563 lines — and carry ``min-w-0`` so a wide
    table inside one cannot drag the row past the page.
    """
    with ui.column().classes("w-full min-w-0 gap-4 p-4"):
        headline = status_headline(cohort)
        with ui.expansion(value=True).classes(
            "w-full min-w-0 border border-gray-200 rounded-md"
        ) as status_expansion:
            with status_expansion.add_slot("header"):
                with ui.row().classes("items-center gap-2 no-wrap w-full min-w-0"):
                    ui.label("Pipeline status").classes(
                        "text-sm font-semibold text-gray-700 whitespace-nowrap"
                    )
                    ui.badge(
                        headline.get("label", ""), color=headline.get("color", "grey")
                    ).classes("text-[10px]")
                    ui.label(headline.get("caption", "")).classes(
                        "text-xs text-gray-500 truncate min-w-0"
                    )
            # The expansion header above already carries the title, the
            # status chip and the step summary.
            render_status_panel(cohort, title=None)

        # flex-wrap so the two panels stack rather than crushing each other on a
        # narrow window; min-w-0 so each can actually shrink to its share.
        with ui.row().classes("w-full min-w-0 gap-4 flex-wrap items-start"):
            with ui.column().classes("flex-1 min-w-0 basis-[24rem]"):
                render_params_panel(cohort)
            with ui.column().classes("flex-1 min-w-0 basis-[24rem]"):
                render_params_doc_panel()


@ui.page("/cohort/{cohort_name}")
def cohort_page(cohort_name: str) -> None:
    """Render the cohort detail page."""
    if redirect := check_auth():
        return redirect
    create_header(cohort_name)

    try:
        store = get_data_store()
        cohort = store.get_cohort(cohort_name)

        if cohort is None:
            # The directory may well exist and simply have no usable pedigree,
            # in which case saying "not found" would send someone looking for
            # the wrong thing.
            stub = next(
                (s for s in store.incomplete_cohorts if s.name == cohort_name), None
            )
            with ui.column().classes("w-full max-w-6xl mx-auto p-6 gap-2"):
                if stub is None:
                    ui.label(f"Cohort not found: {cohort_name}").classes(
                        "text-xl text-red-500"
                    )
                else:
                    ui.label(cohort_name).classes("text-xl font-bold text-amber-900")
                    ui.label(
                        "Unreadable pedigree file"
                        if stub.reason == "unreadable"
                        else "Missing pedigree file"
                    ).classes("text-amber-800 font-semibold")
                    ui.label(
                        "This cohort cannot be explored until"
                        f" {stub.expected_pedigree.name} is in place."
                    ).classes("text-sm text-gray-600")
                    if stub.error:
                        ui.label(stub.error).classes(
                            "text-sm text-amber-900 italic font-mono"
                        )
                ui.button("← Back to Home", on_click=lambda: ui.navigate.to("/"))
            return

        # Load diagnostics per sample (highest priority wins)
        sample_diag: Dict[str, str] = defaultdict(str)
        for fname in ["snvs.tsv", "svs.tsv"]:
            diag_file = store.data_dir / "diagnostics" / fname
            if not diag_file.exists():
                continue
            with open(diag_file, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    if row.get("Ignore", "0") == "1":
                        continue
                    sid = row.get("Sample", "")
                    diag = row.get("Diagnostic", "")
                    if not sid or not diag:
                        continue
                    current = sample_diag[sid]
                    if _DIAG_PRIORITY.get(diag, 0) > _DIAG_PRIORITY.get(current, 0):
                        sample_diag[sid] = diag

        # Build flat list of all individuals from pedigree
        all_individuals: List[Dict[str, Any]] = []
        for family in cohort.families.values():
            for sample in family.samples:
                all_individuals.append(
                    {
                        "FID": sample.family_id,
                        "Sample ID": sample.sample_id,
                        "Sex": sample.sex or "-",
                        "Phenotype": sample.phenotype or "-",
                        "Father": sample.father_id or "-",
                        "Mother": sample.mother_id or "-",
                        "Diagnostic": sample_diag.get(sample.sample_id, ""),
                        "_diag_color": {
                            "pathogenic": "#dc2626",
                            "uncertain": "#f59e0b",
                            "benign": "#16a34a",
                            "conflicting": "#f59e0b",
                        }.get(sample_diag.get(sample.sample_id, ""), ""),
                    }
                )

        # Collect unique values for multiselect filters
        phenotype_values = sorted(
            {ind["Phenotype"] for ind in all_individuals if ind["Phenotype"] != "-"}
        )
        sex_values = sorted(
            {ind["Sex"] for ind in all_individuals if ind["Sex"] != "-"}
        )
        # Diagnostic filter: include actual values + "NA" for undiagnosed
        diag_values_raw = sorted(
            {ind["Diagnostic"] for ind in all_individuals if ind["Diagnostic"]}
        )
        diag_filter_options = diag_values_raw + ["NA"]

        with ui.column().classes("w-full px-6 py-6"):
            # Cohort header
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label(f"🧬 {cohort_name}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{cohort.num_families} families").props("color=blue")
                ui.badge(f"{cohort.num_samples} samples").props("color=teal")

            # Shared filtered state for the right-hand panels
            filtered_state: Dict[str, Any] = {
                "individuals": list(all_individuals),
            }

            # Refresh callbacks of the panels that follow the table filters.
            # Only tabs that have actually been loaded register here, so a
            # filter change never forces work for a tab never opened.
            refreshers: List[Callable[[], None]] = []

            # Page-level tabs, full width. "General" is the individuals table
            # and the panels that follow its filters; "Workflow" is everything
            # about the pipeline run, which belongs to the cohort rather than to
            # the current filter and needs more room than a 400px column.
            has_params = cohort.params_file is not None

            with ui.tabs().classes("w-full") as page_tabs:
                general_tab = ui.tab("General")
                # Absent rather than disabled: with no parameters file there is
                # no pipeline to describe.
                workflow_tab = ui.tab("Workflow") if has_params else None

            # Reading the run record hashes two files and fetching the reference
            # is a network call, so the tab waits until it is opened.
            workflow_state = {"loaded": False}

            with ui.tab_panels(page_tabs, value=general_tab).classes("w-full"):
                with ui.tab_panel(general_tab).classes("w-full p-0"):
                    # Side-by-side layout: table + statistics panel
                    with ui.row().classes("w-full items-start gap-4"):
                        # Left: pedigree table
                        with ui.column():
                            # Count label (updated by filter callback)
                            count_label = ui.label(
                                f"{len(all_individuals)} individuals"
                            ).classes("text-lg font-semibold text-blue-700 mb-2")

                            # Table holder for update_data access
                            dt_ref: Dict[str, Any] = {"dt": None}

                            def on_filter(e: Dict[str, Any]) -> None:
                                filters = e.get("filters", {})
                                filtered = all_individuals

                                # FID text filter
                                fid_text = (filters.get("FID") or "").strip().lower()
                                if fid_text:
                                    filtered = [
                                        ind
                                        for ind in filtered
                                        if fid_text in ind["FID"].lower()
                                    ]

                                # Sample ID text filter
                                sample_text = (
                                    (filters.get("Sample ID") or "").strip().lower()
                                )
                                if sample_text:
                                    filtered = [
                                        ind
                                        for ind in filtered
                                        if sample_text in ind["Sample ID"].lower()
                                    ]

                                # Sex filter (multiselect — list of selected values)
                                sex_vals = filters.get("Sex") or []
                                if sex_vals:
                                    selected_sex = set(sex_vals)
                                    filtered = [
                                        ind
                                        for ind in filtered
                                        if ind["Sex"] in selected_sex
                                    ]

                                # Phenotype filter (multiselect)
                                pheno_vals = filters.get("Phenotype") or []
                                if pheno_vals:
                                    selected = set(pheno_vals)
                                    filtered = [
                                        ind
                                        for ind in filtered
                                        if ind["Phenotype"] in selected
                                    ]

                                # Diagnostic filter (multiselect with NA support)
                                diag_vals = filters.get("Diagnostic") or []
                                if diag_vals:
                                    selected_diag = set(diag_vals)
                                    has_na = "NA" in selected_diag
                                    actual_diag = selected_diag - {"NA"}
                                    filtered = [
                                        ind
                                        for ind in filtered
                                        if (ind["Diagnostic"] in actual_diag)
                                        or (has_na and not ind["Diagnostic"])
                                    ]

                                # Has father checkbox
                                if filters.get("Father"):
                                    filtered = [
                                        ind for ind in filtered if ind["Father"] != "-"
                                    ]

                                # Has mother checkbox
                                if filters.get("Mother"):
                                    filtered = [
                                        ind for ind in filtered if ind["Mother"] != "-"
                                    ]

                                # Update shared filtered state
                                filtered_state["individuals"] = filtered

                                # Update count label
                                label = f"{len(filtered)} individuals"
                                if len(filtered) < len(all_individuals):
                                    label += f" (of {len(all_individuals)} total)"
                                count_label.text = label

                                if dt_ref["dt"]:
                                    dt_ref["dt"].update_data(filtered)

                                # Re-render the panels that follow the filters
                                for refresher in refreshers:
                                    refresher()

                            dt = DataTable(
                                columns=[
                                    {
                                        "id": "FID",
                                        "header": "Family ID",
                                        "cellType": "link",
                                        "href": f"/cohort/{cohort_name}/family/{{FID}}",
                                        "sortable": True,
                                        "minWidth": 250,
                                        "filter": {
                                            "type": "text",
                                            "placeholder": "Filter...",
                                        },
                                    },
                                    {
                                        "id": "Sample ID",
                                        "header": "Sample ID",
                                        "sortable": True,
                                        "minWidth": 90,
                                        "filter": {
                                            "type": "text",
                                            "placeholder": "Filter...",
                                        },
                                    },
                                    {
                                        "id": "Sex",
                                        "header": "Sex",
                                        "sortable": True,
                                        "filter": {
                                            "type": "multiselect",
                                            "options": sex_values,
                                            "placeholder": "All",
                                        },
                                    },
                                    {
                                        "id": "Phenotype",
                                        "header": "Phenotype",
                                        "sortable": True,
                                        "filter": {
                                            "type": "multiselect",
                                            "options": phenotype_values,
                                            "placeholder": "All",
                                        },
                                    },
                                    {
                                        "id": "Father",
                                        "header": "Father ID",
                                        "sortable": True,
                                        "minWidth": 90,
                                        "filter": {
                                            "type": "checkbox",
                                            "label": "Has father",
                                        },
                                    },
                                    {
                                        "id": "Mother",
                                        "header": "Mother ID",
                                        "sortable": True,
                                        "minWidth": 90,
                                        "filter": {
                                            "type": "checkbox",
                                            "label": "Has mother",
                                        },
                                    },
                                    {
                                        "id": "Diagnostic",
                                        "header": "Diagnostic",
                                        "sortable": True,
                                        "cellType": "badge",
                                        "colorField": "_diag_color",
                                        "filter": {
                                            "type": "multiselect",
                                            "options": diag_filter_options,
                                            "placeholder": "All",
                                        },
                                    },
                                ],
                                rows=all_individuals,
                                row_key="Sample ID",
                                pagination={"rowsPerPage": 20},
                                on_filter=on_filter,
                            )
                            dt_ref["dt"] = dt

                        # Right: tabbed panel — Statistics, Ancestry, Polygenic
                        # Scores. Everything describing the pipeline run lives on the
                        # page-level Workflow tab instead: it belongs to the cohort, not
                        # to the current table filters, and needs the full width.
                        has_ancestry = probe_cohort_ancestry_data(cohort)
                        has_pgs = probe_cohort_pgs_data(cohort)

                        # Each tab reads its files on first view, not on first paint.
                        tab_state: Dict[str, Dict[str, bool]] = {
                            "Ancestry": {"loaded": False},
                            "Polygenic Scores": {"loaded": False},
                        }

                        # "min-width: 0" lets the flex item shrink below its content;
                        # without it the tab row and plots push the table sideways.
                        with (
                            ui.column()
                            .classes("flex-1 min-w-[400px]")
                            .style("min-width: 0")
                        ):
                            with ui.tabs().classes("w-full") as panel_tabs:
                                stats_tab = ui.tab("Statistics")
                                ancestry_tab = ui.tab("Ancestry")
                                pgs_tab = ui.tab("Polygenic Scores")
                                if not has_ancestry:
                                    ancestry_tab.props("disable")
                                if not has_pgs:
                                    pgs_tab.props("disable")

                            with ui.tab_panels(panel_tabs, value=stats_tab).classes(
                                "w-full"
                            ):
                                # Each renderer opens its own card, so the panels carry
                                # no border of their own and the chrome is not doubled.
                                with ui.tab_panel(stats_tab).classes("w-full p-0"):
                                    render_stats_panel(store, cohort, filtered_state)

                                with ui.tab_panel(ancestry_tab).classes("w-full p-0"):

                                    @ui.refreshable
                                    def ancestry_content() -> None:
                                        if not has_ancestry:
                                            ui.label(
                                                "No cohort-level ancestry results"
                                            ).classes("text-gray-500 italic")
                                        elif tab_state["Ancestry"]["loaded"]:
                                            render_cohort_ancestry_tab(
                                                store,
                                                cohort,
                                                filtered_state,
                                                refreshers,
                                            )
                                        else:
                                            _loading("Loading ancestry...")

                                    ancestry_content()

                                with ui.tab_panel(pgs_tab).classes("w-full p-0"):

                                    @ui.refreshable
                                    def pgs_content() -> None:
                                        if not has_pgs:
                                            ui.label(
                                                "No cohort-level polygenic score results"
                                            ).classes("text-gray-500 italic")
                                        elif tab_state["Polygenic Scores"]["loaded"]:
                                            render_cohort_pgs_tab(
                                                store,
                                                cohort,
                                                filtered_state,
                                                refreshers,
                                            )
                                        else:
                                            _loading("Loading polygenic scores...")

                                    pgs_content()

                            # Keyed by the exact tab label, and parallel to tab_state:
                            # on_panel_tab_change looks the label up in both, so a tab
                            # in one and not the other is a KeyError or a dead tab.
                            lazy_tabs = {
                                "Ancestry": (has_ancestry, ancestry_content),
                                "Polygenic Scores": (has_pgs, pgs_content),
                            }

                            def on_panel_tab_change(e: Any) -> None:
                                entry = lazy_tabs.get(e.args)
                                if entry is None:
                                    return
                                available, content = entry
                                state = tab_state[e.args]
                                if not available or state["loaded"]:
                                    return

                                def load() -> None:
                                    state["loaded"] = True
                                    content.refresh()

                                ui.timer(0.1, load, once=True)

                            panel_tabs.on("update:model-value", on_panel_tab_change)

                if workflow_tab is not None:
                    with ui.tab_panel(workflow_tab).classes("w-full p-0"):

                        @ui.refreshable
                        def workflow_content() -> None:
                            if not workflow_state["loaded"]:
                                _loading("Loading workflow...")
                                return
                            _render_workflow_tab(cohort)

                        workflow_content()

            def on_page_tab_change(e: Any) -> None:
                if e.args != "Workflow" or workflow_state["loaded"]:
                    return

                async def load() -> None:
                    # Warm both caches off-thread, so the synchronous renderers
                    # only read from cache: one hashes the pedigree and params
                    # file, the other fetches 31 KB over the network, and this
                    # event loop is shared by every session.
                    await asyncio.to_thread(read_cohort_state, cohort.path)
                    await asyncio.to_thread(fetch_params_doc)
                    workflow_state["loaded"] = True
                    workflow_content.refresh()

                ui.timer(0.1, load, once=True)

            page_tabs.on("update:model-value", on_page_tab_change)

    except RuntimeError as e:
        ui.label(f"Error: {e}").classes("text-red-500")
