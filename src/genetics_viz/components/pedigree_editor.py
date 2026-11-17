"""The pedigree editor: a dialog from the family's Sample QC tab, for curators.

Two trees side by side (the pedigree now, and after the changes, the changed
members outlined), the family's rows to edit, the changes worded as ghfc-ngs's
run summary with the rerun they need, the problems found, a required comment
and Save. The history lists earlier edits, each with its trees before and
after and a way to restore the rows before it (saved as a new edit).

The model and the writer are :mod:`genetics_viz.utils.pedigree_edit`; the
plan is the vault note ``Pedigree editor``.
"""

import asyncio
import logging
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from nicegui import ui

from genetics_viz.components.pedigree_tree import (
    OUTLINED_SEVERITIES,
    tree_html,
    worst_severity,
)
from genetics_viz.utils.auth import get_current_user
from genetics_viz.utils.pedigree_edit import (
    ERROR,
    FATHER,
    MOTHER,
    RERUN_CLEAN,
    RERUN_NONE,
    UNKNOWN,
    EditRecord,
    FamilySource,
    Member,
    PedigreeDraft,
    PedigreeEditError,
    describe_changes,
    draft_from_text,
    edit_rows,
    family_cohort_files,
    load_family,
    pending_steps,
    problems,
    read_family_edits,
    rerun_needed,
    save_pedigree_edit_async,
)
from genetics_viz.utils.pedigree_labels import normalize_pheno, normalize_sex
from genetics_viz.utils.pedigree_qc import Fit, QcContext, Suggestion
from genetics_viz.utils.pedigree_svg import PedigreePerson, pedigree_svg
from genetics_viz.utils.sample_qc import FlagRow

logger = logging.getLogger(__name__)

SEX_OPTIONS = {"1": "Male", "2": "Female", UNKNOWN: "Unknown"}
PHENOTYPE_OPTIONS = {"2": "Affected", "1": "Unaffected", UNKNOWN: "Unknown"}
#: The select value that adds a dummy parent for this member.
NEW_DUMMY = "__new_dummy__"

#: What the next ghfc-ngs run has to do, worded for the curator.
RERUN_TEXT = {
    RERUN_NONE: "No ghfc-ngs step reads the phenotype: nothing to rerun.",
    "normal": (
        "Rerun ghfc-ngs as usual: the steps that read the pedigree (Sample QC,"
        " SNVs/Indels and de novo, UPD, STR, genome QC) rebuild this family."
    ),
    RERUN_CLEAN: (
        "The family's sequenced members change: ghfc-ngs needs --clean-stale-families."
    ),
}

#: The changed members' outline in the "after" tree.
CHANGED_COLOR = "#ea580c"
_CHANGED_STYLE = (
    "<style>.gv-pedigree .node.changed .shape"
    f" {{ stroke: {CHANGED_COLOR} !important; stroke-width: 0.12; }}</style>"
)
TREE_PX_PER_UNIT = 60
TREE_MAX_HEIGHT = "320px"


def _parent(code: str) -> str:
    return "" if code == UNKNOWN else code


def changed_ids(original: PedigreeDraft, draft: PedigreeDraft) -> set:
    """Members whose row the edit adds or changes (codes normalised)."""
    was = {m.iid: m.key() for m in [*original.members, *original.implied]}
    return {
        m.iid for m in [*draft.members, *draft.implied] if was.get(m.iid) != m.key()
    }


def draft_people(
    draft: PedigreeDraft,
    flags: Optional[Dict[str, List[FlagRow]]] = None,
    changed: Sequence[str] = (),
) -> List[PedigreePerson]:
    """The draft's members (and implied dummies) as drawing people."""
    people = []
    for m in [*draft.members, *draft.implied]:
        classes: List[str] = []
        if m.iid in changed:
            classes.append("changed")
        elif not m.is_dummy:
            severity = worst_severity((flags or {}).get(m.iid, []))
            if severity in OUTLINED_SEVERITIES:
                classes.append(f"qc-{severity}")
        people.append(
            PedigreePerson(
                m.iid,
                _parent(m.parent(FATHER)),
                _parent(m.parent(MOTHER)),
                normalize_sex(m.sex),
                normalize_pheno(m.phenotype),
                dummy=m.is_dummy,
                classes=tuple(classes),
            )
        )
    return people


def draft_tree_html(
    draft: PedigreeDraft,
    id_prefix: str,
    flags: Optional[Dict[str, List[FlagRow]]] = None,
    changed: Sequence[str] = (),
) -> str:
    """The draft drawn, or a line saying it cannot be."""
    drawing = pedigree_svg(
        draft.family_id, draft_people(draft, flags, changed), id_prefix=id_prefix
    )
    if drawing is None:
        return (
            '<div class="text-sm text-gray-500 italic p-2">'
            "This pedigree cannot be drawn.</div>"
        )
    return _CHANGED_STYLE + tree_html(drawing, TREE_PX_PER_UNIT, TREE_MAX_HEIGHT)


def parent_options(draft: PedigreeDraft, iid: str, role: str) -> Dict[str, str]:
    """A member's father or mother choices: unknown, the family, a new dummy."""
    options = {UNKNOWN: "Unknown"}
    for m in [*draft.members, *draft.implied]:
        if m.iid != iid:
            options[m.iid] = m.iid
    current = draft.get(iid).parent(role) if iid in draft else UNKNOWN
    options.setdefault(current, current)  # a parent with no row
    options[NEW_DUMMY] = f"New dummy {role}"
    return options


def _load(data_dir, family_id: str, files) -> Tuple[FamilySource, List[EditRecord]]:
    return load_family(data_dir, family_id, files), read_family_edits(
        data_dir, family_id
    )


async def open_pedigree_editor(
    anchor: ui.element,
    store: Any,
    family_id: str,
    *,
    flags: Optional[Dict[str, List[FlagRow]]] = None,
    qc: Optional[QcContext] = None,
    on_saved: Optional[Callable[[EditRecord], Any]] = None,
) -> None:
    """Read the family's pedigree off the loop, then open the editor in ``anchor``.

    ``flags`` are the open Sample QC flags by member, for the trees' outlines.
    ``qc`` holds the family's pairs, samples and flags: the draft's fit and the
    suggestions.
    ``on_saved`` gets the edit's record once the files are written.
    """
    files = family_cohort_files(store, family_id)
    try:
        source, history = await asyncio.to_thread(
            _load, store.data_dir, family_id, files
        )
    except PedigreeEditError as e:
        ui.notify(str(e), type="negative", multi_line=True)
        return
    except Exception:
        logger.exception("Failed to read the pedigree of %s", family_id)
        ui.notify(f"The pedigree of {family_id} could not be read", type="negative")
        return
    anchor.clear()  # an earlier editor's dialog
    with anchor:
        _editor_dialog(source, history, flags or {}, qc or QcContext(), on_saved).open()


def _editor_dialog(
    source: FamilySource,
    history: List[EditRecord],
    flags: Dict[str, List[FlagRow]],
    qc: QcContext,
    on_saved: Optional[Callable[[EditRecord], Any]],
) -> ui.dialog:
    family_id = source.family_id
    original = source.original
    state: Dict[str, Any] = {"draft": source.draft(), "undo": [], "answered": set()}
    has_qc = bool(qc.pairs or qc.samples)

    def draft() -> PedigreeDraft:
        return state["draft"]

    def apply(change: Callable[[PedigreeDraft], Any]) -> None:
        before = draft().copy()
        try:
            change(draft())
        except PedigreeEditError as e:
            state["draft"] = before
            ui.notify(str(e), type="warning")
        else:
            state["undo"].append(before)
        refresh()

    with (
        ui.dialog() as dialog,
        ui.card()
        .classes("p-4 gap-3")
        .style("width: min(1400px, 96vw); max-width: 96vw")
        .mark("pedigree-editor"),
    ):
        with ui.row().classes("w-full items-center gap-2"):
            ui.label(f"Edit the pedigree of family {family_id}").classes(
                "text-lg font-semibold text-blue-800"
            )
            ui.label(
                "Writes "
                + ", ".join(f"{c}.pedigree.tsv" for c in source.cohorts)
                + f" and the family's {family_id}.pedigree.tsv"
            ).classes("text-xs text-gray-500")
            ui.space()
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")

        @ui.refreshable
        def trees() -> None:
            changed = changed_ids(original, draft())
            with ui.row().classes("w-full no-wrap gap-4 items-start"):
                for title, which, prefix, marker in (
                    ("Now", original, "ped-edit-now-", "ped-edit-now"),
                    ("After your changes", draft(), "ped-edit-new-", "ped-edit-new"),
                ):
                    with ui.column().classes(
                        "flex-1 min-w-0 border rounded p-2 bg-gray-50 gap-1"
                    ):
                        ui.label(title).classes("text-sm font-semibold text-gray-700")
                        ui.html(
                            draft_tree_html(
                                which,
                                prefix,
                                flags,
                                changed if which is draft() else (),
                            ),
                            sanitize=False,
                        ).classes("w-full").mark(marker)
            ui.label(
                "Orange outline: changed or added. Dashed: not sequenced."
                " Coloured outline: open Sample QC issue."
            ).classes("text-xs text-gray-500")

        trees()

        @ui.refreshable
        def suggestions() -> None:
            if not has_qc:
                return
            found = qc.suggest(draft())
            with (
                ui.column()
                .classes("w-full gap-1 border rounded p-2 bg-blue-50")
                .mark("ped-suggestions")
            ):
                ui.label("Suggestions from Sample QC").classes(
                    "text-sm font-semibold text-blue-900"
                )
                if not found:
                    ui.label(
                        "None: Sample QC suggests no other change for this pedigree."
                    ).classes("text-sm text-gray-600 italic")
                for s in found:
                    _suggestion_row(s, lambda s=s: take(s))

        suggestions()

        with ui.row().classes("w-full items-center gap-2"):
            ui.label("Rows").classes("text-sm font-semibold text-gray-700")
            ui.space()
            ui.button("Undo", icon="undo", on_click=lambda: undo()).props(
                "flat dense size=sm"
            ).mark("ped-undo")
            ui.button("Reset", icon="restart_alt", on_click=lambda: reset()).props(
                "flat dense size=sm"
            ).mark("ped-reset")

        @ui.refreshable
        def rows() -> None:
            d = draft()
            changed = changed_ids(original, d)
            with ui.grid(
                columns="minmax(9rem, 1.4fr) repeat(4, minmax(8rem, 1fr)) auto"
            ).classes("w-full items-center gap-x-3 gap-y-1"):
                for header in ("Member", "Father", "Mother", "Sex", "Phenotype", ""):
                    ui.label(header).classes("text-xs font-semibold text-gray-500")
                for m in [*d.members, *d.implied]:
                    _member_row(d, m, m.iid in changed, m not in d.members, apply)
            _shared_dummy_controls(d, apply)

        rows()

        @ui.refreshable
        def summary() -> None:
            d = draft()
            changes = describe_changes(original, d)
            found = problems(original, d, source.other_ids)
            with ui.row().classes("w-full no-wrap gap-4 items-start"):
                with ui.column().classes("flex-1 min-w-0 gap-1").mark("ped-changes"):
                    ui.label("Changes").classes("text-sm font-semibold text-gray-700")
                    if not changes:
                        ui.label("None yet.").classes("text-sm text-gray-500 italic")
                    for clause in changes:
                        ui.label(clause).classes("text-sm font-mono")
                    if changes:
                        ui.label(RERUN_TEXT[rerun_needed(original, d)]).classes(
                            "text-xs text-blue-800"
                        )
                    if qc.pairs:
                        ui.label(
                            fit_text(qc.fit(original), qc.fit(d), bool(changes))
                        ).classes("text-xs text-gray-700").mark("ped-fit")
                with ui.column().classes("flex-1 min-w-0 gap-1").mark("ped-problems"):
                    ui.label("Problems").classes("text-sm font-semibold text-gray-700")
                    if not found:
                        ui.label("None.").classes("text-sm text-gray-500 italic")
                    for p in found:
                        error = p.level == ERROR
                        with ui.row().classes("no-wrap items-start gap-1"):
                            ui.icon("error" if error else "warning").classes(
                                "text-red-700" if error else "text-amber-700"
                            )
                            ui.label(p.message).classes(
                                "text-sm "
                                + ("text-red-800" if error else "text-amber-900")
                            )
            state["blocked"] = not changes or any(p.level == ERROR for p in found)
            if "update_save" in state:  # not on the first draw: no button yet
                state["update_save"]()

        summary()

        comment = (
            ui.textarea(
                "Why is the pedigree edited? (required)",
                on_change=lambda: update_save(),
            )
            .props("outlined dense autogrow")
            .classes("w-full")
            .mark("ped-comment")
        )
        with ui.row().classes("w-full items-center gap-2"):
            ui.label(
                "The edit is recorded with your name, the date, this comment and"
                " the rows before and after."
            ).classes("text-xs text-gray-500")
            ui.space()
            ui.button("Cancel", on_click=dialog.close).props("flat")
            save_button = (
                ui.button("Save", icon="save", on_click=lambda: save())
                .props("unelevated color=primary")
                .mark("ped-save")
            )

        if history:
            _history(source, history, flags, lambda text: restore(text))

    def update_save() -> None:
        ready = not state.get("blocked", True) and bool((comment.value or "").strip())
        save_button.set_enabled(ready)

    state["update_save"] = update_save

    def refresh() -> None:
        trees.refresh()
        suggestions.refresh()
        rows.refresh()
        summary.refresh()

    def take(s: Suggestion) -> None:
        state["answered"].update(s.flags)
        apply(s.apply)

    def undo() -> None:
        if state["undo"]:
            state["draft"] = state["undo"].pop()
            refresh()

    def reset() -> None:
        if draft().members != original.members:
            state["undo"].append(draft().copy())
            state["draft"] = source.draft()
            state["answered"].clear()
            refresh()

    def restore(text: str) -> None:
        def change(d: PedigreeDraft) -> None:
            restored = draft_from_text(text, family_id)
            d.members = restored.members

        apply(change)
        ui.notify("Rows restored in the draft: review them, then save with a comment")

    async def save() -> None:
        save_button.disable()
        try:
            record = await save_pedigree_edit_async(
                source,
                draft(),
                user=get_current_user(),
                comment=comment.value or "",
                flags=sorted(state["answered"]),
            )
        except (PermissionError, PedigreeEditError) as e:
            ui.notify(str(e), type="negative", multi_line=True)
            update_save()
            return
        except Exception:
            logger.exception("Failed to save the pedigree of %s", family_id)
            ui.notify("The pedigree could not be saved", type="negative")
            update_save()
            return
        ui.notify(
            f"Pedigree of {family_id} saved. " + RERUN_TEXT[record.rerun],
            type="positive",
            multi_line=True,
        )
        dialog.close()
        if on_saved is not None:
            result = on_saved(record)
            if asyncio.iscoroutine(result):
                await result

    update_save()
    return dialog


def _member_row(
    d: PedigreeDraft,
    m: Member,
    changed: bool,
    implied: bool,
    apply: Callable[[Callable[[PedigreeDraft], Any]], None],
) -> None:
    iid = m.iid
    with ui.row().classes("items-center gap-1 no-wrap min-w-0"):
        ui.label(iid).classes(
            "font-mono text-sm truncate "
            + ("text-orange-700 font-semibold" if changed else "")
        )
        if m.is_dummy:
            ui.badge("dummy", color="grey-6").props("outline")
        if implied:
            ui.label("no row").classes("text-xs text-gray-400 italic")

    for role in (FATHER, MOTHER):

        def on_parent(e: Any, role: str = role) -> None:
            value = e.value
            if value == NEW_DUMMY:

                def change(dr: PedigreeDraft) -> None:
                    dr.set_parent(iid, role, None)
                    dr.add_dummy_parents([iid], roles=[role])

                apply(change)
            elif value != m.parent(role):
                apply(lambda dr: dr.set_parent(iid, role, value))

        ui.select(
            parent_options(d, iid, role),
            value=m.parent(role),
            on_change=on_parent,
        ).props("dense outlined options-dense").classes("min-w-0").mark(
            f"ped-{role}-{iid}"
        )

    def on_sex(e: Any) -> None:
        if e.value != normalize_code(m.sex):
            apply(lambda dr: dr.set_sex(iid, e.value))

    def on_phenotype(e: Any) -> None:
        if e.value != normalize_code(m.phenotype):
            apply(lambda dr: dr.set_phenotype(iid, e.value))

    ui.select(SEX_OPTIONS, value=normalize_code(m.sex), on_change=on_sex).props(
        "dense outlined options-dense"
    ).classes("min-w-0").mark(f"ped-sex-{iid}")
    ui.select(
        PHENOTYPE_OPTIONS, value=normalize_code(m.phenotype), on_change=on_phenotype
    ).props("dense outlined options-dense").classes("min-w-0").mark(f"ped-pheno-{iid}")

    with ui.row().classes("items-center gap-0 no-wrap"):
        if not m.is_dummy:
            others = [o for o in d.sequenced_ids if o != iid]
            with (
                ui.button(icon="swap_horiz")
                .props("flat dense size=sm")
                .tooltip("Swap with another member: for a sample swap")
                .mark(f"ped-swap-{iid}")
            ):
                with ui.menu():
                    for other in others:
                        ui.menu_item(
                            f"Swap {iid} and {other}",
                            on_click=lambda o=other: apply(lambda dr: dr.swap(iid, o)),
                        ).mark(f"ped-swap-{iid}-{other}")
        elif not implied:
            ui.button(
                icon="delete", on_click=lambda: apply(lambda dr: dr.remove_dummy(iid))
            ).props("flat dense size=sm").tooltip(
                "Remove this dummy (it must be no one's parent)"
            ).mark(f"ped-remove-{iid}")


def normalize_code(code: str) -> str:
    """``1``, ``2`` or ``0``, as the selects hold them."""
    code = code.strip().removesuffix(".0")
    return code if code in ("1", "2") else UNKNOWN


def _shared_dummy_controls(
    d: PedigreeDraft, apply: Callable[[Callable[[PedigreeDraft], Any]], None]
) -> None:
    """Give several children one dummy parent: full sibs share their parents."""
    lacking = [
        m.iid
        for m in d.members
        if m.parent(FATHER) == UNKNOWN or m.parent(MOTHER) == UNKNOWN
    ]
    if len(lacking) < 2:
        return
    with ui.row().classes("w-full items-center gap-2 pt-2"):
        ui.label("Full sibs:").classes("text-sm text-gray-700")
        children = (
            ui.select(lacking, multiple=True, label="Children")
            .props("dense outlined use-chips options-dense")
            .classes("min-w-[16rem]")
            .mark("ped-shared-children")
        )
        role = (
            ui.select(
                {"both": "Father and mother", FATHER: "Father", MOTHER: "Mother"},
                value="both",
            )
            .props("dense outlined options-dense")
            .mark("ped-shared-role")
        )

        def add() -> None:
            chosen = list(children.value or [])
            if len(chosen) < 2:
                ui.notify("Choose two children or more", type="warning")
                return
            roles = [FATHER, MOTHER] if role.value == "both" else [role.value]
            apply(lambda dr: dr.add_dummy_parents(chosen, roles=roles))

        ui.button("Add shared dummy parents", icon="group_add", on_click=add).props(
            "outline dense"
        ).mark("ped-shared-add")
        ui.label(
            "For the parents they lack. Half sibs: pick one child at a time."
        ).classes("text-xs text-gray-500")


def _history(
    source: FamilySource,
    history: List[EditRecord],
    flags: Dict[str, List[FlagRow]],
    restore: Callable[[str], None],
) -> None:
    with (
        ui.expansion(
            f"History ({len(history)} edit{'s' if len(history) > 1 else ''})",
            icon="history",
        )
        .classes("w-full border rounded")
        .mark("ped-history")
    ):
        for record in reversed(history):
            with ui.column().classes("w-full gap-1 border-b py-2"):
                holder: Dict[str, ui.column] = {}
                with ui.row().classes("w-full items-center gap-2"):
                    ui.label(record.timestamp.replace("T", " ")).classes(
                        "text-sm font-mono"
                    )
                    ui.label(record.user).classes("text-sm font-semibold")
                    ui.label(record.comment).classes("text-sm text-gray-700")
                    ui.space()
                    ui.button(
                        "Show",
                        icon="visibility",
                        on_click=lambda r=record, h=holder: _show_edit(
                            source, r, h["box"], flags, restore
                        ),
                    ).props("flat dense size=sm").mark(f"ped-show-{record.edit_id}")
                for clause in record.changes:
                    ui.label(clause).classes("text-xs font-mono text-gray-600")
                ui.label(RERUN_TEXT.get(record.rerun, record.rerun)).classes(
                    "text-xs text-gray-500"
                )
                holder["box"] = ui.column().classes("w-full")
                holder["box"].set_visibility(False)


async def _show_edit(
    source: FamilySource,
    record: EditRecord,
    box: ui.column,
    flags: Dict[str, List[FlagRow]],
    restore: Callable[[str], None],
) -> None:
    """Toggle an edit's trees before and after, with its Restore button."""
    if box.visible:
        box.set_visibility(False)
        return
    family_id = source.family_id

    def read() -> Tuple[str, str]:
        return (
            edit_rows(source.data_dir, family_id, record.edit_id, "before"),
            edit_rows(source.data_dir, family_id, record.edit_id, "after"),
        )

    try:
        before, after = await asyncio.to_thread(read)
    except OSError:
        ui.notify("The rows kept for this edit are missing", type="negative")
        return
    b = draft_from_text(before, family_id)
    a = draft_from_text(after, family_id)
    changed = changed_ids(b, a)
    box.clear()
    with box:
        with ui.row().classes("w-full no-wrap gap-4 items-start"):
            for title, which, prefix in (
                ("Before", b, f"ped-hist-{record.edit_id}-b-"),
                ("After", a, f"ped-hist-{record.edit_id}-a-"),
            ):
                with ui.column().classes("flex-1 min-w-0 gap-1"):
                    ui.label(title).classes("text-xs text-gray-600")
                    ui.html(
                        draft_tree_html(
                            which, prefix, flags, changed if which is a else ()
                        ),
                        sanitize=False,
                    ).classes("w-full")
        ui.button(
            "Restore the rows before this edit",
            icon="settings_backup_restore",
            on_click=lambda: restore(before),
        ).props("outline dense size=sm").mark(f"ped-restore-{record.edit_id}")
    box.set_visibility(True)


#: ghfc-ngs step names (``#step`` of an ``inputs.tsv``), worded for the tabs.
STEP_NAMES = {
    "sampleqc": "Sample QC",
    "genomeqc": "genome QC",
    "wombat": "SNVs/Indels (Wombat)",
    "denovo": "de novo",
    "upd": "UPD",
    "expansionhunter": "STR",
}


def _pedigree_status(
    store: Any, family_id: str
) -> Tuple[Optional[EditRecord], List[str]]:
    """Blocking: the family's last edit, and the steps still on an older pedigree."""
    files = family_cohort_files(store, family_id)
    edits = read_family_edits(store.data_dir, family_id)
    steps: List[str] = []
    if files:
        text = files[0].read_text(encoding="utf-8")
        found = {p.step for p in pending_steps(store.data_dir, family_id, text)}
        # In the tabs' order, then any step this module does not know.
        steps = [STEP_NAMES[s] for s in STEP_NAMES if s in found]
        steps += sorted(found - set(STEP_NAMES))
    return (edits[-1] if edits else None), steps


def render_pedigree_status(store: Any, family_id: str) -> None:
    """Curators: the last pedigree edit, and the results still on an older pedigree.

    Read off the event loop once the tab is drawn (every ``inputs.tsv`` of the
    family directory). Draws nothing when there is nothing to say.
    """
    box = ui.column().classes("w-full gap-1 mb-2")

    async def load() -> None:
        try:
            last, steps = await asyncio.to_thread(_pedigree_status, store, family_id)
        except Exception:
            logger.warning(
                "Failed to read the pedigree status of %s", family_id, exc_info=True
            )
            return
        box.clear()
        with box:
            if steps:
                with (
                    ui.row()
                    .classes(
                        "w-full items-start no-wrap gap-2 p-2 rounded bg-amber-50"
                        " border border-amber-300"
                    )
                    .mark("pedigree-pending")
                ):
                    ui.icon("update").classes("text-amber-700")
                    ui.label(
                        "Computed on an older pedigree: "
                        + ", ".join(steps)
                        + ". These results change at the next ghfc-ngs run."
                    ).classes("text-sm text-amber-900")
            if last is not None:
                date = last.timestamp.replace("T", " ")[:16]
                label = (
                    ui.label(
                        f"Pedigree last edited on {date} by {last.user}: {last.comment}"
                    )
                    .classes("text-xs text-gray-600")
                    .mark("pedigree-last-edit")
                )
                if last.changes:
                    label.tooltip("; ".join(last.changes))

    ui.timer(0.05, load, once=True)


def _plural(n: int, one: str, many: str) -> str:
    return f"{n} {one if n == 1 else many}"


def fit_text(now: Optional[Fit], after: Optional[Fit], changed: bool) -> str:
    """The Sample QC pairs against the pedigree now and with the changes."""
    if now is None:
        return ""

    def words(f: Fit) -> str:
        text = _plural(f.contradictions, "contradiction", "contradictions")
        if f.unexplained:
            text += ", " + _plural(
                f.unexplained, "pair related beyond it", "pairs related beyond it"
            )
        return text

    text = f"Sample QC against this pedigree: {words(now)} now"
    if changed:
        text += "; " + (words(after) if after else "cannot tell") + " with your changes"
    return text + "."


def _suggestion_row(s: Suggestion, take: Callable[[], None]) -> None:
    with ui.row().classes("w-full items-start no-wrap gap-2"):
        with ui.column().classes("flex-1 min-w-0 gap-0"):
            with ui.row().classes("items-center gap-2"):
                ui.label(s.title).classes("text-sm font-semibold")
                if s.before is not None and s.after is not None:
                    if s.after.contradictions < s.before.contradictions:
                        effect = (
                            f"contradictions {s.before.contradictions}"
                            f" → {s.after.contradictions}"
                        )
                    elif s.after.unexplained < s.before.unexplained:
                        effect = "joins relatives the pedigree does not"
                    else:
                        effect = "settles the sibship"
                    ui.badge(effect, color="blue-1").props("text-color=blue-10")
            for line in s.why:
                ui.label(line).classes("text-xs text-gray-600")
            if s.caution:
                ui.label(s.caution).classes("text-xs text-amber-800")
        ui.button("Apply", on_click=take).props("outline dense size=sm").mark(
            f"ped-apply-{s.key}"
        )
