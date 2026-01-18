"""Home page - displays available cohorts."""

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.utils.data import get_data_store


@ui.page("/")
def home_page() -> None:
    """Render the home/welcome page."""
    create_header()

    with ui.column().classes("w-full max-w-6xl mx-auto p-6"):
        ui.label("Welcome to Genetics-Viz").classes(
            "text-3xl font-bold mb-2 text-blue-900"
        )
        ui.label("Select a cohort to explore").classes("text-lg text-gray-600 mb-6")

        try:
            store = get_data_store()

            if not store.cohorts:
                with ui.card().classes("w-full p-6 bg-yellow-50"):
                    ui.label("⚠️ No cohorts found").classes(
                        "text-xl font-semibold text-yellow-800"
                    )
                    ui.label(
                        f"No valid cohorts were found in: {store.cohorts_dir}"
                    ).classes("text-gray-600")
                    ui.label(
                        "Make sure each cohort directory contains a .pedigree.tsv file."
                    ).classes("text-gray-500 text-sm")
                return

            # Display cohorts as cards
            ui.label("Available Cohorts").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            with ui.row().classes("w-full flex-wrap gap-4"):
                for cohort in sorted(store.cohorts.values(), key=lambda c: c.name):
                    cohort_name = cohort.name  # Capture name for lambda
                    with (
                        ui.card()
                        .classes(
                            "cursor-pointer hover:shadow-lg transition-shadow w-80 border-l-4 border-blue-500"
                        )
                        .on(
                            "click",
                            lambda _, n=cohort_name: ui.navigate.to(f"/cohort/{n}"),
                        )
                    ):
                        with ui.card_section():
                            ui.label(cohort.name).classes(
                                "text-xl font-bold text-blue-700"
                            )

                        with ui.card_section():
                            with ui.row().classes("gap-6"):
                                with ui.column().classes("items-center"):
                                    ui.label(str(cohort.num_families)).classes(
                                        "text-3xl font-bold text-blue-600"
                                    )
                                    ui.label("Families").classes(
                                        "text-sm text-gray-500"
                                    )

                                with ui.column().classes("items-center"):
                                    ui.label(str(cohort.num_samples)).classes(
                                        "text-3xl font-bold text-green-600"
                                    )
                                    ui.label("Samples").classes("text-sm text-gray-500")

                        with ui.card_section().classes("bg-gray-50"):
                            ui.label(f"📄 {cohort.pedigree_file.name}").classes(
                                "text-xs text-gray-400"
                            )

        except RuntimeError as e:
            ui.label(f"Error: {e}").classes("text-red-500")
