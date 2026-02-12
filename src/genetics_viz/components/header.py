"""Header component for genetics-viz."""

from nicegui import ui

from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.score_colors import reload_score_configs
from genetics_viz.utils.gene_scoring import reload_gene_scoring
from genetics_viz.pages.cohort.components.wombat_tab import reload_wombat_configs


def reload_all_configs() -> None:
    """Reload all YAML configuration files."""
    try:
        reload_score_configs()
        reload_gene_scoring()
        reload_wombat_configs()
        ui.notify("✅ Configuration files reloaded successfully", type="positive")
    except Exception as e:
        ui.notify(f"❌ Error reloading configs: {e}", type="negative")


def create_header() -> None:
    """Create the application header with navigation menu."""
    with ui.header().classes("bg-blue-700 text-white items-center justify-between"):
        with ui.row().classes("items-center gap-4"):
            ui.label("🧬 Genetics-Viz").classes("text-xl font-bold")

            # Navigation menu
            with ui.row().classes("gap-2"):
                ui.button(
                    "Home", on_click=lambda: ui.navigate.to("/"), icon="home"
                ).props("flat color=white")

                # Search dropdown menu
                with ui.button("Search", icon="search").props("flat color=white"):
                    with ui.menu():
                        try:
                            store = get_data_store()
                            for cohort_name in sorted(store.cohorts.keys()):
                                ui.menu_item(
                                    cohort_name,
                                    on_click=lambda n=cohort_name: ui.navigate.to(
                                        f"/search/{n}"
                                    ),
                                )
                        except RuntimeError:
                            ui.menu_item("Loading...", auto_close=False)

                # Cohorts dropdown menu
                with ui.button("Cohorts", icon="folder").props("flat color=white"):
                    with ui.menu():
                        try:
                            store = get_data_store()
                            for cohort_name in sorted(store.cohorts.keys()):
                                ui.menu_item(
                                    cohort_name,
                                    on_click=lambda n=cohort_name: ui.navigate.to(
                                        f"/cohort/{n}"
                                    ),
                                )
                        except RuntimeError:
                            ui.menu_item("Loading...", auto_close=False)

                # Validation dropdown menu
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
                            # Separator and bottom items
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

        # Right side - data directory indicator and refresh button
        with ui.row().classes("items-center gap-2"):
            try:
                store = get_data_store()
                ui.label(f"📁 {store.data_dir.name}").classes("text-sm opacity-75")
            except RuntimeError:
                pass

            # Refresh configs button
            ui.button(
                icon="refresh",
                on_click=reload_all_configs
            ).props("flat color=white size=sm round").tooltip("Reload configuration files")
