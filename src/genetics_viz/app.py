"""
NiceGUI web application for genetics-viz.

This is the main entry point for the application. All pages are modular:
- genetics_viz.pages.cohort: Home, Cohort, Family pages
- genetics_viz.pages.validation: Validation Statistics, File, All pages
- genetics_viz.pages.diagnostic: Diagnostic All, Statistics pages
- genetics_viz.pages.search: Cohort-wide search
- genetics_viz.pages.login: Authentication
- genetics_viz.pages.profile: User profile
- genetics_viz.pages.admin: Admin management pages
"""

import os
from pathlib import Path

from nicegui import app as nicegui_app
from nicegui import ui

from genetics_viz.config_model import load_config
from genetics_viz.utils.data import get_static_prefix_map, init_all_data_stores

# Import pages to register routes — this triggers all @ui.page decorators.
from genetics_viz.pages import (  # noqa: F401
    admin,
    cohort,
    diagnostic,
    login,
    profile,
    search,
    validation,
)


def run_app(
    config_file: Path,
    host: str = "127.0.0.1",
    port: int = 8080,
    reload: bool = False,
) -> None:
    """Initialize and run the NiceGUI application.

    Args:
        config_file: Path to the YAML configuration file.
        host: Host address to bind the server to.
        port: Port to run the server on.
        reload: Enable auto-reload for development.
    """
    # Store config in environment for reload mode
    if reload:
        os.environ["GENETICS_VIZ_CONFIG_FILE"] = str(config_file)
        os.environ["GENETICS_VIZ_HOST"] = host
        os.environ["GENETICS_VIZ_PORT"] = str(port)
        os.environ["GENETICS_VIZ_RELOAD"] = "1"

    _init_from_config(config_file)

    config = load_config(config_file)
    ui.run(
        host=host,
        port=port,
        reload=reload,
        title="Genetics-Viz",
        favicon="🧬",
        dark=False,
        storage_secret=config.storage_secret,
    )


def _init_from_config(config_file: Path) -> None:
    """Load configuration and initialize all data stores."""
    config = load_config(config_file)
    init_all_data_stores(config.data_directories)

    # Register static file routes for each data directory
    for path_str, prefix in get_static_prefix_map().items():
        nicegui_app.add_static_files(prefix, path_str)


# Auto-initialize when module is reloaded (for --reload mode)
if os.environ.get("GENETICS_VIZ_RELOAD") == "1":
    config_file_env = os.environ.get("GENETICS_VIZ_CONFIG_FILE")
    if config_file_env:
        _init_from_config(Path(config_file_env))

        config = load_config(Path(config_file_env))
        host = os.environ.get("GENETICS_VIZ_HOST", "127.0.0.1")
        port = int(os.environ.get("GENETICS_VIZ_PORT", "8080"))
        ui.run(
            host=host,
            port=port,
            reload=True,
            title="Genetics-Viz",
            favicon="🧬",
            dark=False,
            storage_secret=config.storage_secret,
        )
