"""Filter components for genetics-viz."""

from typing import Any, Callable, Dict, List

from nicegui import ui

from genetics_viz.components.icons import get_validation_icon


def create_validation_filter_menu(
    all_statuses: List[str],
    filter_state: Dict[str, List[str]],
    on_change: Callable[[], Any],
) -> Dict[str, Any]:
    """Create a validation filter button with dropdown menu.

    Args:
        all_statuses: List of all possible validation statuses
        filter_state: Dictionary with "value" key containing selected statuses
        on_change: Callback function to call when filter changes

    Returns:
        Dictionary mapping status names to their checkbox elements
    """
    validation_checkboxes: Dict[str, Any] = {}

    with (
        ui.button(
            "Filter by Validation",
            icon="filter_list",
        )
        .props("outline color=blue")
        .classes("mb-4")
    ):
        with ui.menu():
            ui.label("Select Validation Statuses:").classes(
                "px-4 py-2 font-semibold text-sm"
            )
            ui.separator()

            with ui.column().classes("p-2"):
                # Quick select buttons
                with ui.row().classes("gap-2 mb-2 flex-wrap"):

                    def select_all():
                        filter_state["value"] = list(all_statuses)
                        for cb in validation_checkboxes.values():
                            cb.value = True
                        on_change()

                    def select_none():
                        filter_state["value"] = []
                        for cb in validation_checkboxes.values():
                            cb.value = False
                        on_change()

                    ui.button(
                        "All",
                        on_click=select_all,
                    ).props("size=sm flat dense").classes("text-xs")
                    ui.button(
                        "None",
                        on_click=select_none,
                    ).props("size=sm flat dense").classes("text-xs")

                ui.separator()

                # Checkboxes for each validation status
                for status in all_statuses:

                    def handle_change(st: str, val: bool) -> None:
                        if val:
                            if st not in filter_state["value"]:
                                filter_state["value"].append(st)
                        else:
                            if st in filter_state["value"]:
                                filter_state["value"].remove(st)
                        on_change()

                    icon, color = get_validation_icon(status)
                    with ui.row().classes("items-center gap-1"):
                        if icon:
                            ui.icon(icon).props(f"size=sm color={color}")
                        validation_checkboxes[status] = ui.checkbox(
                            status,
                            value=status in filter_state["value"],
                            on_change=lambda e, s=status: handle_change(s, e.value),
                        ).classes("text-sm")

    return validation_checkboxes
