"""
Command-line interface for genetics-viz.
"""

from pathlib import Path
from typing import Annotated

import typer

from genetics_viz.app import run_app

app = typer.Typer(
    name="genetics-viz",
    help="A web-based visualization tool for genetics cohort data.",
    add_completion=False,
)


@app.command()
def main(
    config_file: Annotated[
        Path,
        typer.Argument(
            help="Path to the YAML configuration file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            resolve_path=True,
        ),
    ],
    host: Annotated[
        str,
        typer.Option(
            "--host",
            "-h",
            help="Host address to bind the server to",
        ),
    ] = "127.0.0.1",
    port: Annotated[
        int,
        typer.Option(
            "--port",
            "-p",
            help="Port to run the server on",
        ),
    ] = 8080,
    reload: Annotated[
        bool,
        typer.Option(
            "--reload",
            "-r",
            help="Enable auto-reload for development",
        ),
    ] = False,
) -> None:
    """
    Start the genetics-viz web application.

    CONFIG_FILE should be a YAML file listing data directories and users.
    See the README for the expected format.
    """
    typer.echo(f"Starting genetics-viz with config: {config_file}")
    typer.echo(f"Server running at http://{host}:{port}")
    run_app(config_file=config_file, host=host, port=port, reload=reload)


if __name__ in {"__main__", "__mp_main__"}:
    app()
