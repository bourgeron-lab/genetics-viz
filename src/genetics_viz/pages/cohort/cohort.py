"""Cohort detail page - displays families in a cohort."""

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.data import get_data_store


@ui.page("/cohort/{cohort_name}")
def cohort_page(cohort_name: str) -> None:
    """Render the cohort detail page."""
    create_header()

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

        with ui.column().classes("w-full px-6 py-6"):
            # Breadcrumb navigation
            with ui.row().classes("items-center gap-2 mb-4"):
                ui.link("Home", "/").classes("text-blue-600 hover:text-blue-800")
                ui.label("/").classes("text-gray-400")
                ui.label(cohort_name).classes("font-semibold")

            # Cohort header
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label(f"🧬 {cohort_name}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{cohort.num_families} families").props("color=blue")
                ui.badge(f"{cohort.num_samples} samples").props("color=teal")

            # Main content: two-panel layout
            with ui.row().classes("w-full gap-6 flex flex-row"):
                # Left panel: Families table
                with ui.column().classes("flex-1 min-w-0"):
                    ui.label("Families").classes(
                        "text-xl font-semibold mb-2 text-blue-800"
                    )

                    families_data = cohort.get_families_summary()

                    DataTable(
                        columns=[
                            {
                                "id": "Family ID",
                                "header": "Family ID",
                                "cellType": "link",
                                "href": f"/cohort/{cohort_name}/family/{{Family ID}}",
                                "sortable": True,
                            },
                            {
                                "id": "Members",
                                "header": "Members",
                                "sortable": True,
                            },
                        ],
                        rows=families_data,
                        row_key="Family ID",
                        selection="single",
                        pagination={"rowsPerPage": 10},
                        on_selection=lambda e: on_family_select(e),
                    )

                # Right panel: Family members (shown when family selected)
                with ui.column().classes("flex-1 min-w-0"):
                    members_label = ui.label("Select a family to view members").classes(
                        "text-xl font-semibold mb-2 text-gray-400"
                    )
                    members_container = ui.column().classes("w-full")

                def on_family_select(e) -> None:
                    """Handle family selection."""
                    members_container.clear()

                    selected_keys = e.get("selected", [])
                    if not selected_keys:
                        members_label.text = "Select a family to view members"
                        members_label.classes(
                            remove="text-blue-800", add="text-gray-400"
                        )
                        return

                    family_id = e.get("row", {}).get("Family ID")
                    if not family_id:
                        return

                    members_label.text = f"Members of Family: {family_id}"
                    members_label.classes(remove="text-gray-400", add="text-blue-800")

                    members_data = cohort.get_family_members(family_id)

                    with members_container:
                        DataTable(
                            columns=[
                                {"id": "Sample ID", "header": "Sample ID", "sortable": True},
                                {"id": "Father", "header": "Father"},
                                {"id": "Mother", "header": "Mother"},
                                {"id": "Sex", "header": "Sex"},
                                {"id": "Phenotype", "header": "Phenotype"},
                            ],
                            rows=members_data,
                            row_key="Sample ID",
                            pagination={"rowsPerPage": 50},
                        )

    except RuntimeError as e:
        ui.label(f"Error: {e}").classes("text-red-500")
