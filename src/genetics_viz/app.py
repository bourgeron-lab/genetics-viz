"""
NiceGUI web application for genetics-viz.

This is the main entry point for the application. All pages are modular:
- genetics_viz.pages.cohort: Home, Cohort, Family pages
- genetics_viz.pages.validation: Validation Statistics, File, All pages
- genetics_viz.pages.diagnostic: Diagnostic All, Statistics pages
- genetics_viz.pages.search: Cohort-wide search
- genetics_viz.pages.gene: Gene and variant pages, across cohorts
- genetics_viz.pages.login: Authentication
- genetics_viz.pages.profile: User profile
- genetics_viz.pages.admin: Admin management pages
"""

import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from nicegui import app as nicegui_app
from nicegui import ui

from genetics_viz.config_model import load_config
from genetics_viz.utils.change_monitor import start_polling, stop_polling
from genetics_viz.utils.data import (
    get_static_prefix_map,
    init_all_data_stores,
    mount_data_files,
)
from genetics_viz.utils.gene_index import start_gene_index, stop_gene_index
from genetics_viz.utils.genetrek import start_refresh, stop_refresh

# Import pages to register routes — this triggers all @ui.page decorators.
from genetics_viz.pages import (  # noqa: F401
    admin,
    cohort,
    diagnostic,
    gene,
    login,
    profile,
    search,
    validation,
)


@nicegui_app.get("/health")
def health() -> dict:
    """Liveness probe and version report, used by the deployment's status job.

    Deliberately unauthenticated so the F5 pool monitor and CI can reach it.
    Reports nothing beyond liveness and the running release.
    """
    try:
        running = version("genetics-viz")
    except PackageNotFoundError:
        running = "unknown"
    return {"status": "ok", "version": running}


def run_app(
    config_file: Path,
    host: str = "127.0.0.1",
    port: int = 8080,
    reload: bool = False,
    show: bool = False,
) -> None:
    """Initialize and run the NiceGUI application.

    Args:
        config_file: Path to the YAML configuration file.
        host: Host address to bind the server to.
        port: Port to run the server on.
        reload: Enable auto-reload for development.
        show: Open a browser tab on start. Off by default: the app runs as a
            server, where there is no browser to open.
    """
    # Store config in environment for reload mode
    if reload:
        os.environ["GENETICS_VIZ_CONFIG_FILE"] = str(config_file)
        os.environ["GENETICS_VIZ_HOST"] = host
        os.environ["GENETICS_VIZ_PORT"] = str(port)
        os.environ["GENETICS_VIZ_RELOAD"] = "1"
        os.environ["GENETICS_VIZ_SHOW"] = "1" if show else "0"

    _init_from_config(config_file)

    config = load_config(config_file)
    ui.run(
        host=host,
        port=port,
        reload=reload,
        show=show,
        title="Genetics-Viz",
        favicon="🧬",
        dark=False,
        storage_secret=config.storage_secret,
    )


def _init_from_config(config_file: Path) -> None:
    """Load configuration and initialize all data stores."""
    config = load_config(config_file)
    init_all_data_stores(config.data_directories)

    # Serve each data directory's files, to logged-in users only
    for path_str, prefix in get_static_prefix_map().items():
        mount_data_files(prefix, path_str)

    nicegui_app.on_startup(lambda: start_polling(interval=config.poll_interval))
    nicegui_app.on_shutdown(stop_polling)
    # GeneTrek snapshot: the disk cache first, then the current release,
    # both off the event loop. Pages render with neutral badges meanwhile.
    nicegui_app.on_startup(start_refresh)
    nicegui_app.on_shutdown(stop_refresh)
    # Cross-cohort gene index: rebuilt off the event loop where a cohort table
    # changed, served from local disk meanwhile.
    nicegui_app.on_startup(start_gene_index)
    nicegui_app.on_shutdown(stop_gene_index)


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
            show=os.environ.get("GENETICS_VIZ_SHOW") == "1",
            title="Genetics-Viz",
            favicon="🧬",
            dark=False,
            storage_secret=config.storage_secret,
        )
