"""
Entry point for running genetics-viz as a module: python -m genetics_viz
"""

from genetics_viz.cli import app

if __name__ in {"__main__", "__mp_main__"}:
    app()
