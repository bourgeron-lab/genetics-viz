"""Header component for genetics-viz."""

from nicegui import ui

from genetics_viz.utils.clinvar import reload_clinvar_config
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.gene_scoring import reload_gene_scoring
from genetics_viz.utils.score_colors import reload_score_configs
from genetics_viz.utils.vep import reload_vep_config
from genetics_viz.utils.view_presets import reload_view_presets


def reload_all_configs() -> None:
    """Reload all YAML configuration files."""
    try:
        reload_score_configs()
        reload_gene_scoring()
        reload_vep_config()
        reload_clinvar_config()
        reload_view_presets()
        ui.notify("Configuration files reloaded successfully", type="positive")
    except Exception as e:
        ui.notify(f"Error reloading configs: {e}", type="negative")


def create_header(cohort_name: str | None = None) -> None:
    """Create the application header with navigation menu.

    Args:
        cohort_name: The currently active cohort/project name, or None
                     if no cohort is selected (e.g. home, validation pages).
    """
    with ui.header().classes("bg-blue-700 text-white items-center justify-between"):
        with ui.row().classes("items-center gap-4"):
            ui.label("🧬 Genetics-Viz").classes("text-xl font-bold")

            with ui.row().classes("gap-2 items-center"):
                # Home button (always visible)
                ui.button(
                    "Home", on_click=lambda: ui.navigate.to("/"), icon="home"
                ).props("flat color=white")

                # Validation dropdown (always visible)
                with ui.button("Validation", icon="verified").props("flat color=white"):
                    with ui.menu():
                        try:
                            store = get_data_store()
                            to_validate_dir = store.data_dir / "to_validate"
                            if to_validate_dir.exists() and to_validate_dir.is_dir():
                                tsv_files = sorted(
                                    [f.stem for f in to_validate_dir.glob("*.tsv")]
                                )
                                for file_name in tsv_files:
                                    ui.menu_item(
                                        file_name,
                                        on_click=lambda fn=file_name: ui.navigate.to(
                                            f"/validation/file/{fn}"
                                        ),
                                    )
                            ui.separator()
                            ui.menu_item(
                                "See All",
                                on_click=lambda: ui.navigate.to("/validation/all"),
                            )
                            ui.menu_item(
                                "Statistics",
                                on_click=lambda: ui.navigate.to(
                                    "/validation/statistics"
                                ),
                            )
                            ui.separator()
                            ui.menu_item(
                                "Waves",
                                on_click=lambda: ui.navigate.to("/validation/waves"),
                            )
                        except RuntimeError:
                            ui.menu_item("Loading...", auto_close=False)

                # Project selector dropdown (always visible)
                try:
                    store = get_data_store()
                    cohort_names = sorted(store.cohorts.keys())
                except RuntimeError:
                    cohort_names = []

                def on_project_change(e) -> None:
                    if e.value:
                        ui.navigate.to(f"/cohort/{e.value}")

                ui.select(
                    options=cohort_names,
                    value=cohort_name,
                    label="Project",
                    on_change=on_project_change,
                ).props("outlined dense dark color=white label-color=white").classes(
                    "w-48"
                )

                # Cohort button (visible only when a project is selected)
                if cohort_name:
                    ui.button(
                        "Cohort",
                        icon="folder",
                        on_click=lambda n=cohort_name: ui.navigate.to(f"/cohort/{n}"),
                    ).props("flat color=white")

                # Search button (visible only when a project is selected)
                if cohort_name:
                    ui.button(
                        "Search",
                        icon="search",
                        on_click=lambda n=cohort_name: ui.navigate.to(f"/search/{n}"),
                    ).props("flat color=white")

        # Right side - data directory indicator and refresh button
        with ui.row().classes("items-center gap-2"):
            try:
                store = get_data_store()
                ui.label(f"📁 {store.data_dir.name}").classes("text-sm opacity-75")
            except RuntimeError:
                pass

            ui.button(icon="refresh", on_click=reload_all_configs).props(
                "flat color=white size=sm round"
            ).tooltip("Reload configuration files")
