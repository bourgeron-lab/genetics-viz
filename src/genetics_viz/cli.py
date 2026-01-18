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
    data_dir: Annotated[
        Path,
        typer.Argument(
            help="Path to the data directory containing cohort data",
            exists=True,
            file_okay=False,
            dir_okay=True,
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

    DATA_DIR should be the path to a directory containing cohort data,
    with subdirectories under 'cohorts/' each containing a .pedigree.tsv file.
    """
    typer.echo(f"Starting genetics-viz with data from: {data_dir}")
    typer.echo(f"Server running at http://{host}:{port}")
    run_app(data_dir=data_dir, host=host, port=port, reload=reload)


if __name__ in {"__main__", "__mp_main__"}:
    app()
