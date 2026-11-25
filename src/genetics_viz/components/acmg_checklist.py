"""The diagnostic dialog's ACMG criteria checklist (:mod:`genetics_viz.utils.acmg`).

Each criterion is a toggle, grouped by direction and by its own strength. A
ticked one shows its strength, which the curator may change (PM2 is proposed
at Supporting). The points and the class they give are shown live and passed
to ``on_change``, which sets the dialog's verdict; the curator may still choose
another, recorded as an override.
"""

from typing import Callable, Dict, List, Optional

from nicegui import ui

from genetics_viz.utils import acmg
from genetics_viz.utils.diagnostic_badges import verdict_color, verdict_label

#: Rows of the checklist: (direction, strength of the criteria in the row).
_ROWS = [
    ("pathogenic", "VeryStrong"),
    ("pathogenic", "Strong"),
    ("pathogenic", "Moderate"),
    ("pathogenic", "Supporting"),
    ("benign", acmg.STAND_ALONE_STRENGTH),
    ("benign", "Strong"),
    ("benign", "Supporting"),
]
_ON = {"pathogenic": "red-8", "benign": "green-8"}


class AcmgChecklist:
    """Criteria ticked by the curator, their points and class."""

    def __init__(self, on_change: Optional[Callable[[str], None]] = None) -> None:
        self.on_change = on_change
        #: code → strength, for the ticked criteria.
        self.state: Dict[str, str] = {}
        self._refresh: Callable[[], None] = lambda: None

    @property
    def applied(self) -> List[acmg.Applied]:
        return acmg.ordered(acmg.apply(c, s) for c, s in self.state.items())

    @property
    def computed(self) -> str:
        return acmg.classify(self.applied)

    def render(self) -> None:
        with (
            ui.expansion("ACMG criteria (optional)", icon="checklist")
            .classes("w-full border rounded")
            .props("dense")
            .mark("acmg-criteria")
        ):
            ui.label(
                "Richards et al. 2015 criteria, added up as points (Tavtigian et "
                "al. 2020): pathogenic from 10, likely pathogenic 6 to 9, "
                "uncertain 0 to 5, likely benign -1 to -6, benign from -7. "
                "Hover a criterion for its definition."
            ).classes("text-xs text-gray-500")
            grid = ui.refreshable(self._grid)
            grid()
            self._refresh = grid.refresh

    def _grid(self) -> None:
        applied = self.applied
        points = acmg.total_points(applied)
        computed = acmg.classify(applied)
        with ui.row().classes("items-center gap-2 mt-1").mark("acmg-points"):
            if applied:
                ui.label(f"{points} points").classes("text-sm font-semibold")
                ui.label("→").classes("text-gray-400")
                ui.label(verdict_label(computed)).classes(
                    "text-sm font-semibold"
                ).style(f"color: {verdict_color(computed)}")
                ui.label(acmg.format_criteria(applied).replace(";", ", ")).classes(
                    "text-xs text-gray-500"
                )
                ui.button("Clear", icon="backspace", on_click=self._clear).props(
                    "flat dense no-caps size=sm"
                ).mark("acmg-clear")
            else:
                ui.label("No criterion applied").classes("text-sm text-gray-500 italic")
        for direction, strength in _ROWS:
            codes = [
                c.code
                for c in acmg.CRITERIA.values()
                if c.direction == direction and c.strength == strength
            ]
            if not codes:
                continue
            with ui.row().classes("items-center gap-1 w-full flex-wrap"):
                ui.label(
                    f"{direction.capitalize()}, {acmg.STRENGTH_LABELS[strength].lower()}"
                ).classes("text-xs text-gray-500 w-36 shrink-0")
                for code in codes:
                    self._toggle(code)

    def _toggle(self, code: str) -> None:
        criterion = acmg.CRITERIA[code]
        on = code in self.state
        with ui.row().classes("items-center gap-0 no-wrap"):
            button = (
                ui.button(code, on_click=lambda c=code: self._flip(c))
                .props("dense no-caps size=sm")
                .mark(f"acmg-{code}")
            )
            if on:
                button.props(f"unelevated color={_ON[criterion.direction]}")
            else:
                button.props("outline color=grey-7")
            button.tooltip(criterion.text)
            if on and len(criterion.strengths) > 1:
                # A menu, not a select: a dense select is taller than the
                # toggles and made the rows uneven.
                color = _ON[criterion.direction]
                with (
                    ui.button(acmg.STRENGTH_LABELS[self.state[code]])
                    .props(
                        f"flat dense no-caps size=sm color={color}"
                        " icon-right=arrow_drop_down"
                    )
                    .classes("ml-0.5 mr-2")
                    .mark(f"acmg-{code}-strength")
                ):
                    with ui.menu():
                        for strength in criterion.strengths:
                            ui.menu_item(
                                acmg.STRENGTH_LABELS[strength],
                                on_click=lambda c=code, st=strength: self._set_strength(
                                    c, st
                                ),
                            ).props("dense").mark(f"acmg-{code}-{strength}")

    def _flip(self, code: str) -> None:
        if code in self.state:
            del self.state[code]
        else:
            self.state[code] = acmg.CRITERIA[code].default
        self._changed()

    def _set_strength(self, code: str, strength: str) -> None:
        if strength and code in self.state and strength != self.state[code]:
            self.state[code] = strength
            self._changed()

    def _clear(self) -> None:
        self.state.clear()
        self._changed()

    def _changed(self) -> None:
        self._refresh()
        if self.on_change:
            self.on_change(self.computed)
