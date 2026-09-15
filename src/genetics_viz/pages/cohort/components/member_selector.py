"""Family member selection table - shared by the cohort and standalone family pages."""

from typing import Any, Callable, Dict, List

from nicegui import ui

# A single CSS grid holds the header and every data cell as direct children, so
# all columns resolve against one shared set of tracks. Applying the grid style
# per row instead makes each row an independent grid whose `auto`/`1fr` tracks
# size to that row's own content only, which leaves the headers misaligned with
# the values below them.
# `gap: 0` keeps the row stripes contiguous - horizontal spacing comes from the
# per-cell `px-2` padding instead. The `max-content` floor on the ID columns
# stops a header such as "Sample ID" wrapping onto two lines when its `1fr`
# share is narrower than the text; the parent scrolls if the table cannot fit.
_GRID_STYLE = (
    "display: grid;"
    " grid-template-columns: auto auto minmax(max-content, 1fr)"
    " minmax(max-content, 1fr) minmax(max-content, 1fr) auto auto;"
    " gap: 0;"
    # Cells stretch to the row height so the stripes and bottom borders line up
    # across the row - without this the empty header above the "only" buttons is
    # only as tall as its padding and leaves a notch in the header strip.
    " align-items: stretch;"
)

# Second column holds the "only" buttons and needs no header.
_HEADERS = ("Select", "", "Sample ID", "Father", "Mother", "Sex", "Phenotype")

_PEDIGREE_FIELDS = ("Father", "Mother", "Sex", "Phenotype")

_CELL = "px-2 py-1"

# A stretched cell is taller than the control it holds, so centre it explicitly.
_CONTROL_CELL_STYLE = "display: flex; align-items: center"


def render_member_selector(
    members_data: List[Dict[str, Any]],
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the family pedigree table with per-member display checkboxes.

    Toggling a member updates ``selected_members["value"]`` in place and runs
    every callback in ``data_table_refreshers``, so the dependent panels and
    variant tables re-filter.

    :param members_data: member dicts as returned by ``Cohort.get_family_members``
    :param selected_members: mutable ``{"value": [sample_id, ...]}`` selection state
    :param data_table_refreshers: refresh callbacks to run on every selection change
    """
    member_checkboxes: Dict[str, ui.checkbox] = {}

    def refresh_all() -> None:
        for refresher in data_table_refreshers:
            refresher()

    def select_all_members() -> None:
        selected_members["value"] = [m["Sample ID"] for m in members_data]
        for checkbox in member_checkboxes.values():
            checkbox.value = True
        refresh_all()

    def select_none_members() -> None:
        selected_members["value"] = []
        for checkbox in member_checkboxes.values():
            checkbox.value = False
        refresh_all()

    def make_change_handler(sid):
        def handler(e):
            if e.value and sid not in selected_members["value"]:
                selected_members["value"].append(sid)
            elif not e.value and sid in selected_members["value"]:
                selected_members["value"].remove(sid)
            refresh_all()

        return handler

    def make_only_handler(sid):
        def handler():
            selected_members["value"] = [sid]
            for s_id, checkbox in member_checkboxes.items():
                checkbox.value = s_id == sid
            refresh_all()

        return handler

    # `tight()` drops the card's own 1rem padding and gap so the table can sit
    # flush inside the tinted panel below.
    with ui.card().tight().classes("flex-1"):
        with (
            ui.column()
            .classes("w-full p-3 bg-blue-50")
            .style("gap: 0.5rem; overflow-x: auto")
        ):
            with ui.row().classes("items-center gap-2"):
                ui.label("Select Members to Display:").classes(
                    "font-semibold text-blue-800"
                )
                ui.button("All", on_click=select_all_members).props(
                    "size=sm flat dense"
                ).classes("text-xs")
                ui.button("None", on_click=select_none_members).props(
                    "size=sm flat dense"
                ).classes("text-xs")

            with ui.element("div").classes("w-full text-sm").style(_GRID_STYLE):
                for header in _HEADERS:
                    ui.label(header).classes(f"{_CELL} bg-blue-100 font-semibold")

                # There is no per-row wrapper element: the seven cells of a row
                # share a background and a bottom border, which reads as a row.
                for idx, member in enumerate(members_data):
                    sample_id = member["Sample ID"]
                    cell = (
                        f"{_CELL} border-b border-gray-200 "
                        f"{'bg-white' if idx % 2 == 0 else 'bg-gray-50'}"
                    )

                    with ui.element("div").classes(cell).style(_CONTROL_CELL_STYLE):
                        member_checkboxes[sample_id] = ui.checkbox(
                            "",
                            value=True,
                            on_change=make_change_handler(sample_id),
                        ).props("dense size=xs")
                    with ui.element("div").classes(cell).style(_CONTROL_CELL_STYLE):
                        ui.button("only", on_click=make_only_handler(sample_id)).props(
                            "size=xs flat dense color=blue"
                        ).classes("text-xs")
                    ui.label(sample_id).classes(f"{cell} font-medium")
                    for field in _PEDIGREE_FIELDS:
                        ui.label(member.get(field, "-")).classes(
                            f"{cell} text-gray-600"
                        )
