"""Filter components for genetics-viz."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Dict, List, Tuple

from nicegui import ui

from genetics_viz.components.icons import get_validation_icon

if TYPE_CHECKING:
    from genetics_viz.utils.genesets import GenesetSource


def create_validation_filter_menu(
    all_statuses: List[str],
    filter_state: Dict[str, List[str]],
    on_change: Callable[[], Any],
    label: str = "Filter by Validation",
    icon_fn: Callable[[str], Tuple[str, str]] | None = None,
    button_classes: str = "mb-4",
    button_size: str = "",
    button_props: str = "",
) -> Dict[str, Any]:
    """Create a validation filter button with dropdown menu.

    Args:
        all_statuses: List of all possible validation statuses
        filter_state: Dictionary with "value" key containing selected statuses
        on_change: Callback function to call when filter changes
        label: Button label text
        button_classes: CSS classes for the button
        button_size: Button size (e.g., "h-14")
        button_props: Extra Quasar props for the button (e.g., "dense")

    Returns:
        Dictionary mapping status names to their checkbox elements
    """
    validation_checkboxes: Dict[str, Any] = {}

    props_str = "outline color=blue"
    if button_props:
        props_str += f" {button_props}"
    button = ui.button(
        label,
        icon="filter_list",
    ).props(props_str)

    if button_classes:
        button.classes(button_classes)
    if button_size:
        button.classes(button_size)

    with button:
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

                    icon, color = (icon_fn or get_validation_icon)(status)
                    with ui.row().classes("items-center gap-1"):
                        if icon:
                            ui.icon(icon).props(f"size=sm color={color}")
                        validation_checkboxes[status] = ui.checkbox(
                            status,
                            value=status in filter_state["value"],
                            on_change=lambda e, s=status: handle_change(s, e.value),
                        ).classes("text-sm")

    return validation_checkboxes


def create_geneset_filter_menu(
    selected: Dict[str, List[str]],
    sources: List[GenesetSource],
    on_change: Callable[[], Any] | None = None,
    button_props: str = "",
) -> None:
    """A "Genesets" button whose menu selects gene sets, grouped by source.

    ``selected["value"]`` holds the selected ``GenesetSource.key``s. The
    button turns green while anything is selected. GeneTrek's older versions
    and sub-lists are behind a switch, unless already in use.
    """
    if not sources:
        return
    known = {s.key for s in sources}
    selected["value"] = [k for k in selected["value"] if k in known]
    view = {"search": "", "more": False}

    button = ui.button("Genesets", icon="list")

    def sync_button() -> None:
        if selected["value"]:
            button.props(remove="outline", add="unelevated color=green")
        else:
            button.props(remove="unelevated color=green", add="outline")
        button.update()

    button.props(button_props)
    sync_button()

    def changed() -> None:
        sync_button()
        if on_change is not None:
            on_change()

    def set_many(keys: List[str], value: bool) -> None:
        current = set(selected["value"])
        current = current | set(keys) if value else current - set(keys)
        # Keep the menu's order, so the selection reads predictably.
        selected["value"] = [s.key for s in sources if s.key in current]
        listing.refresh()
        changed()

    @ui.refreshable
    def listing() -> None:
        needle = view["search"].strip().lower()
        groups: Dict[str, List[Any]] = {}
        for s in sources:
            if needle and needle not in s.label.lower() and needle not in s.key.lower():
                continue
            if s.secondary and not view["more"] and not needle:
                if s.key not in selected["value"]:
                    continue
            groups.setdefault(s.group, []).append(s)
        if not groups:
            ui.label("No matching geneset").classes("text-sm text-gray-500 p-2")
        for group, items in groups.items():
            n_sel = sum(s.key in selected["value"] for s in items)
            title = f"{group} ({n_sel}/{len(items)})" if n_sel else group
            with (
                ui.expansion(title, value=bool(n_sel) or bool(needle))
                .props("dense")
                .classes("w-full")
            ):
                with ui.row().classes("gap-1"):
                    keys = [s.key for s in items]
                    ui.button(
                        "All", on_click=lambda _, k=keys: set_many(k, True)
                    ).props("size=sm flat dense").classes("text-xs")
                    ui.button(
                        "None", on_click=lambda _, k=keys: set_many(k, False)
                    ).props("size=sm flat dense").classes("text-xs")
                for s in items:
                    cb = ui.checkbox(
                        f"{s.label} ({s.n_genes:,} genes)",
                        value=s.key in selected["value"],
                        on_change=lambda e, k=s.key: set_many([k], e.value),
                    ).classes("text-sm")
                    if s.description:
                        cb.tooltip(s.description)

    with button:
        with ui.menu().props("max-height=70vh"):
            with ui.column().classes("p-2 gap-1 w-96"):
                with ui.row().classes("w-full items-center gap-2 no-wrap"):

                    def on_search(e) -> None:
                        view["search"] = e.value or ""
                        listing.refresh()

                    ui.input(placeholder="Search genesets", on_change=on_search).props(
                        "dense outlined clearable"
                    ).classes("grow")
                    ui.button(
                        "Clear", on_click=lambda: set_many(list(known), False)
                    ).props("size=sm flat dense")

                def on_more(e) -> None:
                    view["more"] = e.value
                    listing.refresh()

                if any(s.secondary for s in sources):
                    ui.switch("Older versions and sub-lists", on_change=on_more).props(
                        "dense"
                    ).classes("text-xs")
                listing()
