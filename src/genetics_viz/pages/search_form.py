"""Search page form: the source component cards and the Individuals panel.

The form only writes into its state (:class:`SearchForm`); the page reads it
when Search is pressed. Logic without widgets is in ``utils/search_rows.py``.
"""

import re
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

from nicegui import ui

from genetics_viz.components.filters import (
    create_consequence_filter_menu,
    create_geneset_filter_menu,
    create_validation_filter_menu,
)
from genetics_viz.utils.search_rows import (
    ALL_CALL_VALUES,
    IndividualFilters,
    new_wcx_component,
    new_wombat_component,
    source_key,
)
from genetics_viz.utils.validation_status import FILTER_STATUSES

_LOCUS_PLACEHOLDER = "chr1:10000-10100, SHANK3, ENSG..."


@dataclass
class SearchForm:
    """What the form edits, and what it needs to offer."""

    wombat_sources: List[Dict[str, Any]]
    wcx_sources: List[Dict[str, Any]]
    genesets: List[Any]
    sex_values: List[str]
    phenotype_values: List[str]
    components: List[Dict[str, Any]] = field(default_factory=list)
    individuals: IndividualFilters = field(default_factory=IndividualFilters)
    #: The locus inputs on screen, rebuilt with the component list.
    locus_inputs: List[Any] = field(default_factory=list)
    #: Set once the page defines its search: Enter in a locus box runs it.
    on_search: Optional[Callable[..., Any]] = None

    def new_wombat(self) -> Dict[str, Any]:
        default = source_key(self.wombat_sources[0]) if self.wombat_sources else ""
        return new_wombat_component(default)

    def new_wcx(self) -> Dict[str, Any]:
        default = source_key(self.wcx_sources[0]) if self.wcx_sources else ""
        return new_wcx_component(default)

    def reset_components(self) -> None:
        """One default wombat component, as the page opens."""
        self.components.clear()
        if self.wombat_sources:
            self.components.append(self.new_wombat())

    def bind_locus_enter(self, loc_input: Any) -> None:
        """Bind Enter key on a locus input to trigger search."""
        if self.on_search:
            loc_input.on("keydown.enter", self.on_search)

    def set_search(self, fn: Callable[..., Any]) -> None:
        """Install the search and bind Enter on the locus inputs already built.

        Inputs built later (an added component) are bound as they are built.
        """
        self.on_search = fn
        for loc_inp in self.locus_inputs:
            loc_inp.on("keydown.enter", fn)


def render_search_form(form: SearchForm) -> Any:
    """Build the search panel. Returns the Search button."""
    with ui.card().classes("w-full p-2 mb-2").props("flat bordered"):
        # Side-by-side panels
        with ui.row().classes("w-full gap-4 items-start flex-nowrap"):
            # --- LEFT PANEL: Variants (2/3 width) ---
            with ui.column().classes("flex-[2]"):
                ui.label("Variants").classes("text-sm font-semibold text-gray-600")
                render_component_list = _component_list(form)
                render_component_list()

            ui.separator().props("vertical")

            # --- RIGHT PANEL: Individuals (1/3 width) ---
            with ui.column().classes("flex-[1] min-w-0 w-full"):
                ui.label("Individuals").classes("text-sm font-semibold text-gray-600")
                clear_individuals = _individuals_panel(form)

        # Search + Clear buttons centered below both panels
        with ui.row().classes("justify-center gap-4 mt-2"):
            search_button = (
                ui.button("Search", icon="search")
                .props("color=blue dense")
                .mark("search-button")
            )

            def _clear_all():
                form.reset_components()
                form.individuals.clear()
                clear_individuals()
                render_component_list.refresh()

            ui.button("Clear", icon="clear_all", on_click=_clear_all).props(
                "outline dense"
            )
    return search_button


def _component_list(form: SearchForm):
    """The refreshable list of component cards, with its add buttons."""

    def remove(comp_id: str) -> None:
        form.components[:] = [c for c in form.components if c["id"] != comp_id]
        render_component_list.refresh()

    @ui.refreshable
    def render_component_list():
        form.locus_inputs.clear()
        for i, comp in enumerate(form.components):
            if i > 0:
                ui.label("OR").classes(
                    "text-xs font-bold text-center w-full text-gray-400 my-0 py-0"
                )
            if comp["type"] == "wombat":
                _render_wombat_card(form, comp, remove)
            elif comp["type"] == "wisecondorx":
                _render_wcx_card(form, comp, remove)

        # Add buttons row
        with ui.row().classes("gap-2 mt-1"):
            if form.wombat_sources:

                def _add_wombat():
                    form.components.append(form.new_wombat())
                    render_component_list.refresh()

                ui.button("+ Wombat", on_click=_add_wombat).props(
                    "flat dense size=sm color=blue no-caps"
                )
            if form.wcx_sources:

                def _add_wcx():
                    form.components.append(form.new_wcx())
                    render_component_list.refresh()

                ui.button("+ WisecondorX", on_click=_add_wcx).props(
                    "flat dense size=sm color=teal no-caps"
                )

    return render_component_list


def _locus_input(form: SearchForm, comp: Dict[str, Any]) -> None:
    loc_input = (
        ui.input(
            label="Locus (optional)",
            placeholder=_LOCUS_PLACEHOLDER,
            value=comp["locus"]["value"],
            on_change=lambda e, c=comp: c["locus"].update({"value": e.value or ""}),
        )
        .props("outlined dense")
        .classes("flex-grow")
    )
    form.locus_inputs.append(loc_input)
    form.bind_locus_enter(loc_input)


def _remove_button(comp: Dict[str, Any], remove: Callable[[str], None]) -> None:
    ui.button(icon="close", on_click=lambda _, i=comp["id"]: remove(i)).props(
        "flat round dense size=sm color=grey"
    )


def _box_checkbox(label: str, box: Dict[str, Any]) -> None:
    """A small checkbox writing into a ``{"value": bool}`` box."""
    ui.checkbox(
        label,
        value=box["value"],
        on_change=lambda e, b=box: b.update({"value": e.value}),
    ).props("dense").classes("text-xs")


def _validation_menu(comp: Dict[str, Any]) -> None:
    create_validation_filter_menu(
        all_statuses=list(FILTER_STATUSES),
        filter_state=comp["validations"],
        on_change=lambda: None,
        label="Validation",
        button_classes="",
        button_props="dense size=sm",
    )


def _render_wombat_card(
    form: SearchForm, comp: Dict[str, Any], remove: Callable[[str], None]
) -> None:
    """Render a compact wombat source filter card."""
    with ui.card().classes("w-full p-2 mb-0").props("flat bordered"):
        # Row 1: label + file selector + locus + remove
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("Wombat").classes("text-sm font-semibold text-blue-700")
            ui.select(
                options={
                    source_key(sf): sf["display_name"] for sf in form.wombat_sources
                },
                value=comp["source_key"],
                label="Source",
                on_change=lambda e, c=comp: c.update({"source_key": e.value}),
            ).props("outlined dense").classes("w-64")
            _locus_input(form, comp)
            _remove_button(comp, remove)

        # Row 2: genesets + impacts + validation + exclude checkboxes
        with ui.row().classes("items-center gap-2 w-full flex-wrap"):
            create_geneset_filter_menu(
                comp["genesets"], form.genesets, button_props="dense size=sm"
            )
            create_consequence_filter_menu(
                comp["impacts"], button_props="dense size=sm"
            )
            _validation_menu(comp)
            _box_checkbox("Exclude LCR", comp["exclude_lcr"])
            _box_checkbox("Exclude gnomAD filtered", comp["exclude_gnomad"])
            _box_checkbox("Exclude gnomAD WGS", comp["exclude_gnomad_wgs"])


def _render_wcx_card(
    form: SearchForm, comp: Dict[str, Any], remove: Callable[[str], None]
) -> None:
    """Render a compact wisecondorx source filter card."""
    with ui.card().classes("w-full p-2 mb-0").props("flat bordered"):
        # Row 1: label + exonic only + locus + remove
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("WisecondorX").classes("text-sm font-semibold text-teal-700")
            _box_checkbox("Exonic only", comp["exonic_only"])
            _locus_input(form, comp)
            _remove_button(comp, remove)

        # Row 2: genesets + impacts + validation + ratio filter
        with ui.row().classes("items-center gap-2 w-full flex-wrap"):
            create_geneset_filter_menu(
                comp["genesets"], form.genesets, button_props="dense size=sm"
            )
            _call_filter_menu(comp)
            _validation_menu(comp)

            ui.label("|ratio|:").classes("text-xs text-gray-600")
            ui.input(
                label="Min",
                value=str(comp["ratio_min"]["value"] or ""),
                on_change=_ratio_handler(comp, "ratio_min"),
            ).props("outlined dense").classes("w-20")
            ui.input(
                label="Max",
                value=str(comp["ratio_max"]["value"] or ""),
                on_change=_ratio_handler(comp, "ratio_max"),
            ).props("outlined dense").classes("w-20")


def _ratio_handler(comp: Dict[str, Any], key: str):
    """Parse a ``|ratio|`` bound into ``comp[key]``; a bad number marks the input."""

    def handler(e):
        raw = (e.value or "").replace(",", ".")
        if not raw:
            comp[key].update({"value": None})
            e.sender._props["error"] = False
            e.sender.update()
            return
        try:
            comp[key].update({"value": float(raw)})
            e.sender._props["error"] = False
            e.sender.update()
        except ValueError:
            e.sender._props["error"] = True
            e.sender.update()
            ui.notify(
                f"Invalid number: {raw}",
                type="warning",
                position="top",
                timeout=2000,
            )

    return handler


def _call_filter_menu(comp: Dict[str, Any]) -> None:
    """Build an Impacts dropdown menu for a wisecondorx component."""
    comp_calls = comp["selected_calls"]
    call_btn = ui.button("Impacts", icon="filter_list").props("outline dense size=sm")
    call_cbs: Dict[str, Any] = {}

    def _update_call_btn():
        if len(comp_calls["value"]) == len(ALL_CALL_VALUES):
            call_btn.props(remove="unelevated color=orange", add="outline")
        else:
            call_btn.props(remove="outline", add="unelevated color=orange")
        call_btn.update()

    def _set_all(value: bool):
        comp_calls["value"] = list(ALL_CALL_VALUES) if value else []
        for cb in call_cbs.values():
            cb.value = value
        _update_call_btn()

    def make_call_handler(cv):
        def handler(e):
            if e.value:
                if cv not in comp_calls["value"]:
                    comp_calls["value"].append(cv)
            else:
                if cv in comp_calls["value"]:
                    comp_calls["value"].remove(cv)
            _update_call_btn()

        return handler

    with call_btn:
        with ui.menu():
            ui.label("Filter by Impact:").classes("px-4 py-2 font-semibold text-sm")
            ui.separator()
            with ui.column().classes("p-2"):
                with ui.row().classes("gap-2 mb-2"):
                    ui.button("All", on_click=lambda: _set_all(True)).props(
                        "size=sm flat dense"
                    ).classes("text-xs")
                    ui.button("None", on_click=lambda: _set_all(False)).props(
                        "size=sm flat dense"
                    ).classes("text-xs")

                ui.separator()
                for call_val in ALL_CALL_VALUES:
                    call_cbs[call_val] = ui.checkbox(
                        call_val,
                        value=call_val in comp_calls["value"],
                        on_change=make_call_handler(call_val),
                    ).classes("text-sm")


def _individuals_panel(form: SearchForm) -> Callable[[], None]:
    """Sex, phenotype, both-parents and excluded samples.

    Returns a function that empties the widgets (the Clear button), after
    ``form.individuals.clear()`` has emptied the state.
    """
    filters = form.individuals
    with ui.column().classes("gap-2 w-full flex-1"):
        # Sex + Phenotype on the same row
        with ui.row().classes("w-full gap-2"):
            sex_select = (
                ui.select(
                    options=form.sex_values,
                    label="Sex",
                    value=filters.sex["value"],
                    multiple=True,
                    on_change=lambda e: filters.sex.update({"value": e.value or []}),
                )
                .props("outlined dense use-chips")
                .classes("min-w-0")
                .style("flex: 1;")
            )
            phenotype_select = (
                ui.select(
                    options=form.phenotype_values,
                    label="Phenotype",
                    value=filters.phenotype["value"],
                    multiple=True,
                    on_change=lambda e: filters.phenotype.update(
                        {"value": e.value or []}
                    ),
                )
                .props("outlined dense use-chips")
                .classes("min-w-0")
                .style("flex: 2;")
            )

        has_parents_cb = ui.checkbox(
            "Only samples with both parents",
            value=filters.has_parents["value"],
            on_change=lambda e: filters.has_parents.update({"value": e.value}),
        ).props("dense")

        # Exclude samples — chips + input in a
        # bordered container that stretches vertically
        ui.label("Exclude samples").classes("text-xs text-gray-500 -mb-2")
        with (
            ui.column()
            .classes("w-full flex-1 border rounded p-2 gap-1 min-h-[60px]")
            .style("border-color: #bdbdbd;")
        ):
            # Chips area (cleared & rebuilt, input
            # stays outside so it is never destroyed)
            exclude_chips_row = ui.row().classes("w-full gap-1 flex-wrap")
            excluded = filters.exclude_samples

            def _rebuild_exclude_chips():
                """Rebuild chip elements only."""
                exclude_chips_row.clear()
                if excluded["value"]:
                    with exclude_chips_row:
                        for sid in excluded["value"]:

                            def make_remove(s):
                                def remove(_):
                                    excluded["value"] = [
                                        x for x in excluded["value"] if x != s
                                    ]
                                    _rebuild_exclude_chips()

                                return remove

                            ui.chip(
                                sid,
                                removable=True,
                                on_value_change=make_remove(sid),
                                color="red",
                            ).props("dense outline size=sm")

            def _commit_exclude_text(client_value: str | None = None):
                """Parse input text into chips."""
                raw = (
                    client_value
                    if client_value is not None
                    else (exclude_input.value or "")
                )
                tokens = [t.strip() for t in re.split(r"[,\s]+", raw) if t.strip()]
                if tokens:
                    current = excluded["value"]
                    for t in tokens:
                        if t not in current:
                            current.append(t)
                    excluded["value"] = current
                    exclude_input.value = ""
                    _rebuild_exclude_chips()

            def _on_exclude_event(e):
                val = ""
                if isinstance(e.args, dict):
                    val = e.args.get("value", "")
                _commit_exclude_text(val)

            # Persistent input (never destroyed)
            exclude_input = (
                ui.input(placeholder="Type sample IDs...")
                .props("dense hide-bottom-space")
                .classes("w-full exclude-input-wrap")
            )
            exclude_input.on(
                "keydown",
                _on_exclude_event,
                js_handler="""(e) => {
                    if (e.key === ' ' || e.key === ',' || e.key === 'Enter') {
                        e.preventDefault();
                        const el = e.target;
                        emit({value: el.value || ''});
                    }
                }""",
                throttle=0.1,
            )
            exclude_input.on(
                "blur",
                _on_exclude_event,
                js_handler="""(e) => {
                    const el = e.target;
                    const val = el.value || '';
                    if (val.trim()) {
                        emit({value: val});
                    }
                }""",
            )

    def clear_widgets() -> None:
        sex_select.value = []
        phenotype_select.value = []
        has_parents_cb.value = False
        exclude_input.value = ""
        _rebuild_exclude_chips()

    return clear_widgets
