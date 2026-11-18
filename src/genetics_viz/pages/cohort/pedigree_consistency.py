"""``/pedigree-consistency``: the data directory's pedigree disagreements (curators)."""

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.pedigree_consistency import render_all_findings
from genetics_viz.utils.auth import can_write, check_auth
from genetics_viz.utils.data import get_data_store


@ui.page("/pedigree-consistency")
async def pedigree_consistency_page() -> None:
    if redirect := check_auth():
        return redirect
    create_header()
    with ui.column().classes("w-full px-6 py-6 gap-3"):
        ui.label("Pedigrees across cohorts").classes("text-2xl font-bold text-blue-900")
        if not can_write():
            ui.label("Curators and administrators only.").classes("text-gray-600")
            return
        ui.label(
            "The cohorts of a data directory share one families/ directory, and"
            " ghfc-ngs builds a family from the rows of the cohort it runs. A"
            " family two cohorts list differently is rebuilt by each, and one"
            " sample under two family IDs is one family analysed twice."
        ).classes("text-sm text-gray-700 max-w-4xl")
        render_all_findings(get_data_store())
