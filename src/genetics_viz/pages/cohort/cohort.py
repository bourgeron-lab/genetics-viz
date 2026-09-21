"""Cohort detail page - displays pedigree members with filters."""

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
from genetics_viz.pages.cohort.components.stats_panel import render_stats_panel
from genetics_viz.utils.auth import check_auth
from genetics_viz.utils.data import get_data_store

# Diagnostic priority: higher value wins
_DIAG_PRIORITY = {"pathogenic": 3, "uncertain": 2, "benign": 1}


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
            with ui.column().classes("w-full max-w-6xl mx-auto p-6"):
                ui.label(f"Cohort not found: {cohort_name}").classes(
                    "text-xl text-red-500"
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
                        sample_text = (filters.get("Sample ID") or "").strip().lower()
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
                                ind for ind in filtered if ind["Sex"] in selected_sex
                            ]

                        # Phenotype filter (multiselect)
                        pheno_vals = filters.get("Phenotype") or []
                        if pheno_vals:
                            selected = set(pheno_vals)
                            filtered = [
                                ind for ind in filtered if ind["Phenotype"] in selected
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
                            filtered = [ind for ind in filtered if ind["Father"] != "-"]

                        # Has mother checkbox
                        if filters.get("Mother"):
                            filtered = [ind for ind in filtered if ind["Mother"] != "-"]

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

                # Right: tabbed panel — Statistics, Ancestry, Polygenic Scores
                has_ancestry = probe_cohort_ancestry_data(cohort)
                has_pgs = probe_cohort_pgs_data(cohort)

                # Each tab reads its TSVs on first view, not on first paint.
                tab_state: Dict[str, Dict[str, bool]] = {
                    "Ancestry": {"loaded": False},
                    "Polygenic Scores": {"loaded": False},
                }

                # "min-width: 0" lets the flex item shrink below its content;
                # without it the tab row and plots push the table sideways.
                with ui.column().classes("flex-1 min-w-[400px]").style("min-width: 0"):
                    with ui.tabs().classes("w-full") as panel_tabs:
                        stats_tab = ui.tab("Statistics")
                        ancestry_tab = ui.tab("Ancestry")
                        pgs_tab = ui.tab("Polygenic Scores")
                        if not has_ancestry:
                            ancestry_tab.props("disable")
                        if not has_pgs:
                            pgs_tab.props("disable")

                    with ui.tab_panels(panel_tabs, value=stats_tab).classes("w-full"):
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
                                        store, cohort, filtered_state, refreshers
                                    )
                                else:
                                    with ui.column().classes(
                                        "w-full items-center justify-center py-16"
                                    ):
                                        ui.spinner(size="xl", color="blue")
                                        ui.label("Loading ancestry...").classes(
                                            "text-lg text-gray-600 mt-4"
                                        )

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
                                        store, cohort, filtered_state, refreshers
                                    )
                                else:
                                    with ui.column().classes(
                                        "w-full items-center justify-center py-16"
                                    ):
                                        ui.spinner(size="xl", color="blue")
                                        ui.label("Loading polygenic scores...").classes(
                                            "text-lg text-gray-600 mt-4"
                                        )

                            pgs_content()

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

    except RuntimeError as e:
        ui.label(f"Error: {e}").classes("text-red-500")
