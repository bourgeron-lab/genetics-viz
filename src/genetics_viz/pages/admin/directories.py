"""Admin page for managing data directories."""

from pathlib import Path

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.config_model import (
    DataDirectoryConfig,
    get_config,
    save_config,
)
from genetics_viz.utils.auth import check_auth, is_admin
from genetics_viz.utils.data import (
    add_data_store,
    remove_data_store,
)


@ui.page("/admin/directories")
def admin_directories_page() -> None:
    """Manage data directories."""
    if redirect := check_auth():
        return redirect
    if not is_admin():
        ui.navigate.to("/")
        return
    create_header()

    with ui.column().classes("w-full max-w-4xl mx-auto p-6"):
        ui.label("Manage Data Directories").classes(
            "text-3xl font-bold text-blue-900 mb-6"
        )

        # Directory list
        dir_container = ui.column().classes("w-full mb-6")

        def refresh_dirs() -> None:
            dir_container.clear()
            config = get_config()
            with dir_container:
                if not config.data_directories:
                    ui.label("No data directories configured").classes(
                        "text-gray-500 italic"
                    )
                    return

                for i, d in enumerate(config.data_directories):
                    dir_path = Path(d.path)
                    exists = dir_path.exists()

                    with ui.card().classes("w-full p-4 mb-2"):
                        with ui.row().classes("items-center justify-between w-full"):
                            with ui.column().classes("flex-1"):
                                with ui.row().classes("items-center gap-2"):
                                    ui.label(d.path).classes(
                                        "font-mono text-sm font-bold"
                                    )
                                    if d.default:
                                        ui.badge("default", color="green").classes(
                                            "text-xs"
                                        )
                                    if not exists:
                                        ui.badge("not found", color="red").classes(
                                            "text-xs"
                                        )
                                if d.description:
                                    ui.label(d.description).classes(
                                        "text-sm text-gray-500"
                                    )

                            with ui.row().classes("items-center gap-2"):
                                if not d.default:

                                    def make_set_default(idx):
                                        def handler():
                                            cfg = get_config()
                                            for dd in cfg.data_directories:
                                                dd.default = False
                                            cfg.data_directories[idx].default = True
                                            save_config(cfg)
                                            refresh_dirs()
                                            ui.notify(
                                                "Default directory updated",
                                                type="positive",
                                            )

                                        return handler

                                    ui.button(
                                        "Set Default",
                                        icon="star",
                                        on_click=make_set_default(i),
                                    ).props("flat color=blue size=sm")

                                if len(config.data_directories) > 1:

                                    def make_remove(path_str):
                                        def handler():
                                            cfg = get_config()
                                            cfg.data_directories = [
                                                dd
                                                for dd in cfg.data_directories
                                                if dd.path != path_str
                                            ]
                                            # Ensure there's still a default
                                            if cfg.data_directories and not any(
                                                dd.default
                                                for dd in cfg.data_directories
                                            ):
                                                cfg.data_directories[0].default = True
                                            save_config(cfg)
                                            try:
                                                remove_data_store(path_str)
                                            except KeyError:
                                                pass
                                            refresh_dirs()
                                            ui.notify(
                                                "Directory removed", type="positive"
                                            )

                                        return handler

                                    ui.button(
                                        icon="delete",
                                        on_click=make_remove(d.path),
                                    ).props("flat color=red size=sm round")

        refresh_dirs()

        # Add directory form
        with ui.card().classes("w-full p-6"):
            ui.label("Add Data Directory").classes("text-xl font-semibold mb-4")
            path_input = (
                ui.input("Directory path")
                .props("outlined dense")
                .classes("w-full mb-2")
            )
            desc_input = (
                ui.input("Description (optional)")
                .props("outlined dense")
                .classes("w-full mb-4")
            )

            def add_directory() -> None:
                path_str = path_input.value.strip()
                desc = desc_input.value.strip()

                if not path_str:
                    ui.notify("Please enter a directory path", type="warning")
                    return

                dir_path = Path(path_str)
                if not dir_path.exists():
                    ui.notify(f"Directory does not exist: {path_str}", type="negative")
                    return

                config = get_config()
                if any(d.path == path_str for d in config.data_directories):
                    ui.notify("Directory already configured", type="warning")
                    return

                # Add to config
                is_first = len(config.data_directories) == 0
                config.data_directories.append(
                    DataDirectoryConfig(
                        path=path_str,
                        description=desc,
                        default=is_first,
                    )
                )
                save_config(config)

                # Hot-add data store
                try:
                    add_data_store(path_str, desc)
                    ui.notify("Directory added successfully", type="positive")
                except Exception as e:
                    ui.notify(
                        f"Added to config but failed to load: {e}", type="warning"
                    )

                path_input.value = ""
                desc_input.value = ""
                refresh_dirs()

            ui.button("Add Directory", icon="add", on_click=add_directory).props(
                "color=blue unelevated"
            )
