"""Reusable column selector dialog with grouped layout and presets."""

from collections import OrderedDict
from typing import Any, Callable, Dict, List, Tuple

from nicegui import ui

from genetics_viz.utils.column_names import get_column_group, get_display_label


def build_column_selector(
    all_columns: List[str],
    selected_cols: Dict[str, Any],
    on_visibility_change: Callable[[], None],
    presets: List[Dict[str, Any]],
) -> Tuple[ui.dialog, Callable[[], None]]:
    """Build a column selector dialog with grouped layout and presets.

    Parameters
    ----------
    all_columns : list of column IDs available in the data
    selected_cols : mutable dict ``{"value": [...]}``, modified in-place
    on_visibility_change : called after any selection change
    presets : list of preset dicts with "name" and "columns" keys

    Returns
    -------
    (dialog, sync_fn) – open dialog with ``dialog.open()``;
    call ``sync_fn()`` to re-sync checkbox states from external changes.
    """
    # Build group structure (preserve insertion order from all_columns)
    groups: OrderedDict[str, List[str]] = OrderedDict()
    ungrouped: List[str] = []

    for col in all_columns:
        group = get_column_group(col)
        if group:
            groups.setdefault(group, []).append(col)
        else:
            ungrouped.append(col)

    # Checkbox references
    col_cbs: Dict[str, ui.checkbox] = {}
    group_cbs: Dict[str, ui.checkbox] = {}
    _syncing = {"value": False}

    def _sync_all():
        """Sync every checkbox to match selected_cols."""
        _syncing["value"] = True
        for col, cb in col_cbs.items():
            cb.value = col in selected_cols["value"]
        for gname, cb in group_cbs.items():
            gcols = groups[gname]
            cb.value = all(c in selected_cols["value"] for c in gcols)
        _syncing["value"] = False

    def _reorder_and_apply():
        selected_cols["value"] = [
            c for c in all_columns if c in selected_cols["value"]
        ]
        on_visibility_change()
        _sync_all()

    def _handle_col(col_name: str, is_checked: bool):
        if _syncing["value"]:
            return
        if is_checked and col_name not in selected_cols["value"]:
            selected_cols["value"].append(col_name)
        elif not is_checked and col_name in selected_cols["value"]:
            selected_cols["value"].remove(col_name)
        _reorder_and_apply()

    def _handle_group(group_name: str, is_checked: bool):
        if _syncing["value"]:
            return
        for col in groups[group_name]:
            if is_checked and col not in selected_cols["value"]:
                selected_cols["value"].append(col)
            elif not is_checked and col in selected_cols["value"]:
                selected_cols["value"].remove(col)
        _reorder_and_apply()

    def _select_all():
        selected_cols["value"] = list(all_columns)
        _reorder_and_apply()

    def _select_none():
        selected_cols["value"] = []
        _reorder_and_apply()

    def _make_preset_handler(preset):
        def handler():
            available = [
                c for c in preset.get("columns", []) if c in all_columns
            ]
            selected_cols["value"] = available
            _reorder_and_apply()
        return handler

    # --- Build dialog UI ---
    with ui.dialog() as dialog, ui.card().classes(
        "w-[900px] max-w-[95vw]"
    ).style("padding: 12px 16px;"):
        # Header + All/None + presets — single compact row
        with ui.row().classes("w-full items-center gap-2 flex-wrap"):
            ui.label("Columns").classes("text-sm font-semibold mr-1")
            ui.button("All", on_click=_select_all).props(
                "size=xs outline dense"
            )
            ui.button("None", on_click=_select_none).props(
                "size=xs outline dense"
            )
            if presets:
                ui.separator().props("vertical").classes("h-5")
                for preset in presets:
                    label = preset["name"].replace(" View", "")
                    ui.button(
                        label, on_click=_make_preset_handler(preset)
                    ).props("size=xs flat dense").classes("text-xs")
            ui.space()
            ui.button(icon="close", on_click=dialog.close).props(
                "flat round dense size=xs"
            )

        ui.separator().classes("my-1")

        # Column checkboxes — dense CSS multi-column layout
        with ui.element("div").style(
            "column-count: 4; column-gap: 1rem;"
        ):
            # Ungrouped columns
            for col in ungrouped:
                col_cbs[col] = ui.checkbox(
                    get_display_label(col),
                    value=col in selected_cols["value"],
                    on_change=lambda e, c=col: _handle_col(c, e.value),
                ).props("dense").classes("text-xs w-full")

            # Grouped columns
            for group_name, group_cols in groups.items():
                with ui.element("div").style(
                    "break-inside: avoid;"
                ).classes("mt-2"):
                    all_checked = all(
                        c in selected_cols["value"] for c in group_cols
                    )
                    group_cbs[group_name] = ui.checkbox(
                        group_name,
                        value=all_checked,
                        on_change=lambda e, g=group_name: _handle_group(
                            g, e.value
                        ),
                    ).props("dense").classes("text-xs font-bold")

                    with ui.column().classes("pl-5 gap-0"):
                        for col in group_cols:
                            col_cbs[col] = ui.checkbox(
                                get_display_label(col),
                                value=col in selected_cols["value"],
                                on_change=lambda e, c=col: _handle_col(
                                    c, e.value
                                ),
                            ).props("dense").classes("text-xs")

    return dialog, _sync_all
