"""The member card's Tree view: the family drawn by pypedigree.

The member selector (``pages/cohort/components/member_selector.py``) swaps its
table for this drawing. Selection stays in the table: the tree only shows it,
the members not selected dimmed. Each sequenced member's outline takes the
colour of their worst open Sample QC flag (high or medium), as the table's
warning icon does. Dummies, and a parent with no row, are dashed.

The drawing is a static SVG (:func:`pedigree_svg.pedigree_svg`): ``ui.html``
sets ``innerHTML``, which never runs a script. It is redrawn whole on every
change (about 1 ms), so nothing is looked up by ``data-id`` in JavaScript.
The page's tree and the Enlarge dialog's are two drawings with their own id
prefix.
"""

from typing import Any, Callable, Dict, List, Optional, Sequence

from nicegui import ui

from genetics_viz.utils.pedigree_labels import normalize_pheno, normalize_sex
from genetics_viz.utils.pedigree_svg import (
    PedigreeDrawing,
    PedigreePerson,
    pedigree_svg,
)
from genetics_viz.utils.sample_qc import FlagRow, severity_info

#: On the page: pixels per layout unit (a symbol is 0.8 unit, so 44 px; labels
#: about 11 px), and the height the tree may take before it shrinks to fit.
PX_PER_UNIT = 55
MAX_HEIGHT_PX = 300
#: In the Enlarge dialog: larger, within the window.
BIG_PX_PER_UNIT = 90

#: Severities drawn as an outline colour; ``info`` is not an issue to show.
OUTLINED_SEVERITIES = ("high", "medium")

_MISSING = {"", "-", "0", "-9"}


def _parent(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text in _MISSING else text


def worst_severity(flags: Sequence[FlagRow]) -> Optional[str]:
    """The highest severity among a member's open flags, ``None`` without one."""
    if not flags:
        return None
    return max(flags, key=lambda f: f.rank).severity


def tree_people(
    members_data: List[Dict[str, Any]],
    selected: Sequence[str],
    qc_flags: Optional[Dict[str, List[FlagRow]]] = None,
) -> List[PedigreePerson]:
    """The selector's member dicts as drawing people, classed for the tree's CSS."""
    people = []
    for m in members_data:
        sid = m["Sample ID"]
        dummy = bool(m.get("Dummy"))
        classes: List[str] = []
        if not dummy:
            if sid not in selected:
                classes.append("deselected")
            severity = worst_severity((qc_flags or {}).get(sid, []))
            if severity in OUTLINED_SEVERITIES:
                classes.append(f"qc-{severity}")
        people.append(
            PedigreePerson(
                sid,
                _parent(m.get("Father")),
                _parent(m.get("Mother")),
                normalize_sex(m.get("Sex")),
                normalize_pheno(m.get("Phenotype")),
                dummy=dummy,
                classes=tuple(classes),
            )
        )
    return people


def _style() -> str:
    rules = [
        # pypedigree writes ``background: white`` inline; the card is tinted.
        ".gv-pedigree svg { display: block; width: 100%; height: auto;"
        " background: transparent !important; }",
        ".gv-pedigree .node.deselected { opacity: 0.3; }",
    ]
    for severity in OUTLINED_SEVERITIES:
        color = severity_info(severity)["color"]
        rules.append(
            f".gv-pedigree .node.qc-{severity} .shape"
            f" {{ stroke: {color}; stroke-width: 0.09; }}"
        )
    return "<style>" + " ".join(rules) + "</style>"


def tree_html(drawing: PedigreeDrawing, px_per_unit: float, max_height: str) -> str:
    """The drawing in a box sized from its layout, capped by width and height.

    ``max_height`` is a CSS length (``300px``, ``75vh``); the box's width is
    the smallest of the card's, the natural size and what fits that height.
    """
    svg = drawing.svg.split("?>", 1)[-1].strip()
    ratio = drawing.width / drawing.height if drawing.height else 1
    width = (
        f"min(100%, {drawing.width * px_per_unit:.0f}px,"
        f" calc({max_height} * {ratio:.4f}))"
    )
    return (
        f'{_style()}<div class="gv-pedigree" style="width: {width}; margin: 0 auto">'
        f"{svg}</div>"
    )


def tree_key(people: Sequence[PedigreePerson], drawing: PedigreeDrawing) -> List[str]:
    """The key under the tree: only what this drawing shows."""
    parts = []
    if any("deselected" in p.classes for p in people):
        parts.append("Faded: not selected")
    if drawing.has_unsampled:
        parts.append("Dashed: not sequenced")
    shown = [
        severity_info(s)["label"].lower()
        for s in OUTLINED_SEVERITIES
        if any(f"qc-{s}" in p.classes for p in people)
    ]
    if shown:
        parts.append(f"Coloured outline: open Sample QC issue ({', '.join(shown)})")
    return parts


def render_pedigree_tree(
    members_data: List[Dict[str, Any]],
    selected_members: Dict[str, List[str]],
    qc_flags: Optional[Dict[str, List[FlagRow]]] = None,
) -> Callable[[Optional[Dict[str, List[FlagRow]]]], None]:
    """Draw the tree and its Enlarge dialog; return ``redraw(qc_flags=None)``.

    ``redraw`` follows the current selection; given flags, it keeps them for
    later redraws too (after a Sample QC review).
    """
    state: Dict[str, Any] = {"flags": qc_flags or {}}

    def drawing(id_prefix: str):
        people = tree_people(members_data, selected_members["value"], state["flags"])
        return people, pedigree_svg("family", people, id_prefix=id_prefix)

    with (
        ui.dialog() as dialog,
        ui.card().classes("items-center").style("max-width: 95vw; width: auto"),
    ):
        with ui.row().classes("w-full items-center justify-between"):
            ui.label("Pedigree").classes("font-semibold text-blue-800")
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")
        big = ui.html("", sanitize=False).classes("w-full").mark("member-tree-big")

    def fill_big() -> None:
        _, d = drawing("member-tree-big-")
        if d is not None:
            big.set_content(tree_html(d, BIG_PX_PER_UNIT, "80vh"))
            big.style("width: 90vw")
        dialog.open()

    with ui.column().classes("w-full items-stretch").style("gap: 0.25rem"):
        canvas = ui.html("", sanitize=False).classes("w-full").mark("member-tree")
        with ui.row().classes("w-full items-start justify-between no-wrap gap-2"):
            key = ui.label("").classes("text-xs text-gray-500")
            ui.button(icon="open_in_full", on_click=fill_big).props(
                "flat dense size=sm"
            ).tooltip("Enlarge the pedigree").mark("member-tree-enlarge")

    def redraw(flags: Optional[Dict[str, List[FlagRow]]] = None) -> None:
        if flags is not None:
            state["flags"] = flags
        people, d = drawing("member-tree-")
        if d is None:
            canvas.set_content(
                '<div class="text-sm text-gray-500 italic p-2">'
                "The pedigree could not be drawn.</div>"
            )
            key.set_text("")
            return
        canvas.set_content(tree_html(d, PX_PER_UNIT, f"{MAX_HEIGHT_PX}px"))
        parts = tree_key(people, d)
        key.set_text(". ".join(parts) + "." if parts else "")

    redraw()
    return redraw
