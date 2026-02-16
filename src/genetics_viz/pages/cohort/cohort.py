"""Cohort detail page - displays pedigree members with filters."""

from typing import Any, Dict, List

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.pages.cohort.components.stats_panel import render_stats_panel
from genetics_viz.utils.data import get_data_store


@ui.page("/cohort/{cohort_name}")
def cohort_page(cohort_name: str) -> None:
    """Render the cohort detail page."""
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
                    }
                )

        # Collect unique phenotype values for multiselect options
        phenotype_values = sorted(
            {ind["Phenotype"] for ind in all_individuals if ind["Phenotype"] != "-"}
        )

        # Sex options for multiselect filter
        sex_values = sorted(
            {ind["Sex"] for ind in all_individuals if ind["Sex"] != "-"}
        )

        with ui.column().classes("w-full px-6 py-6"):
            # Cohort header
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label(f"🧬 {cohort_name}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{cohort.num_families} families").props("color=blue")
                ui.badge(f"{cohort.num_samples} samples").props("color=teal")

            # Shared filtered state for statistics panel
            filtered_state: Dict[str, Any] = {
                "individuals": list(all_individuals),
            }

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
                                ind for ind in filtered
                                if fid_text in ind["FID"].lower()
                            ]

                        # Sample ID text filter
                        sample_text = (
                            filters.get("Sample ID") or ""
                        ).strip().lower()
                        if sample_text:
                            filtered = [
                                ind for ind in filtered
                                if sample_text in ind["Sample ID"].lower()
                            ]

                        # Sex filter (multiselect — list of selected values)
                        sex_vals = filters.get("Sex") or []
                        if sex_vals:
                            selected_sex = set(sex_vals)
                            filtered = [
                                ind for ind in filtered
                                if ind["Sex"] in selected_sex
                            ]

                        # Phenotype filter (multiselect)
                        pheno_vals = filters.get("Phenotype") or []
                        if pheno_vals:
                            selected = set(pheno_vals)
                            filtered = [
                                ind for ind in filtered
                                if ind["Phenotype"] in selected
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
                        ],
                        rows=all_individuals,
                        row_key="Sample ID",
                        pagination={"rowsPerPage": 20},
                        on_filter=on_filter,
                    )
                    dt_ref["dt"] = dt

                # Right: statistics panel
                with ui.column().classes("flex-1 min-w-[400px]"):
                    render_stats_panel(store, cohort, filtered_state)

    except RuntimeError as e:
        ui.label(f"Error: {e}").classes("text-red-500")
