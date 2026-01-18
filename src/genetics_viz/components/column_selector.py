"""Shared column selector menu component."""

from typing import Any, Callable, Dict, List

from nicegui import ui


def create_column_selector_menu(
    all_columns: List[str],
    selected_cols: Dict[str, List[str]],
    on_change: Callable[[], None],
    get_display_label: Callable[[str], str] | None = None,
) -> None:
    """Create a column selector menu button.

    Args:
        all_columns: List of all available column names
        selected_cols: Dict with 'value' key containing list of selected columns
        on_change: Callback to call when selection changes
        get_display_label: Optional function to get display label for column
    """
    with ui.button("Select Columns", icon="view_column").props("outline color=blue"):
        with ui.menu():
            ui.label("Show/Hide Columns:").classes("px-4 py-2 font-semibold text-sm")
            ui.separator()

            with ui.column().classes("p-2"):
                # Select All / Deselect All buttons
                with ui.row().classes("gap-2 mb-2"):
                    checkboxes: Dict[str, Any] = {}

                    def select_all():
                        selected_cols["value"] = list(all_columns)
                        for cb in checkboxes.values():
                            cb.value = True
                        on_change()

                    def select_none():
                        selected_cols["value"] = []
                        for cb in checkboxes.values():
                            cb.value = False
                        on_change()

                    ui.button("All", on_click=select_all).props(
                        "size=sm flat dense"
                    ).classes("text-xs")
                    ui.button("None", on_click=select_none).props(
                        "size=sm flat dense"
                    ).classes("text-xs")

                ui.separator()

                # Checkboxes for each column
                for col in all_columns:
                    display_name = get_display_label(col) if get_display_label else col

                    def make_handler(column: str):
                        def handler(e):
                            if e.value and column not in selected_cols["value"]:
                                selected_cols["value"].append(column)
                            elif not e.value and column in selected_cols["value"]:
                                selected_cols["value"].remove(column)
                            on_change()

                        return handler

                    checkboxes[col] = ui.checkbox(
                        display_name,
                        value=col in selected_cols["value"],
                        on_change=make_handler(col),
                    ).classes("text-sm")
