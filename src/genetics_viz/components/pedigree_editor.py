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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Set, Tuple

from nicegui import ui

from genetics_viz.components.pedigree_tree import (
    OUTLINED_SEVERITIES,
    tree_html,
    worst_severity,
)
from genetics_viz.utils.auth import get_current_user
from genetics_viz.utils.cohort_state import read_cohort_state
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
    cross_problems,
    describe_changes,
    draft_from_text,
    edit_rows,
    family_cohort_files,
    load_family,
    move_member,
    outside_ids,
    pending_steps,
    read_edit_index,
    read_family_edits,
    WaitingFamily,
    rerun_needed,
    save_cross_family_edit_async,
    save_pedigree_edit_async,
    swap_across,
    waiting_families,
)
from genetics_viz.utils.pedigree_labels import normalize_pheno, normalize_sex
from genetics_viz.utils.pedigree_merge import (
    MergeSource,
    Orphan,
    orphaned_families,
    curation_counts,
    load_merge,
    merge_changes,
    merge_problems,
    merge_rerun,
    save_merge_async,
)
from genetics_viz.utils.pedigree_qc import (
    CrossSuggestion,
    Fit,
    QcContext,
    Suggestion,
    combined_draft,
    merge_candidates,
    merged_context,
    suggest_across,
)
from genetics_viz.utils.pedigree_svg import PedigreePerson, pedigree_svg
from genetics_viz.utils.sample_qc_acks import load_acks
from genetics_viz.utils.sample_qc import (
    FlagRow,
    PairRow,
    flags_by_sample,
    load_cohort_qc,
    load_family_qc,
    pair_rows,
)

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


@dataclass
class Side:
    """One family of the editor: what was read, the draft, its Sample QC."""

    source: FamilySource
    draft: PedigreeDraft
    qc: QcContext
    flags: Dict[str, List[FlagRow]]

    @property
    def family_id(self) -> str:
        return self.source.family_id

    @property
    def original(self) -> PedigreeDraft:
        return self.source.original


def _other_family(store: Any, family_id: str) -> Tuple[FamilySource, QcContext]:
    """Blocking: a second family to edit beside the first, with its own QC."""
    source = load_family(
        store.data_dir, family_id, family_cohort_files(store, family_id)
    )
    data = load_family_qc(store.data_dir, family_id)
    if data is None:
        return source, QcContext()
    return source, QcContext(
        pairs=[p for p in pair_rows(data.pairs) if not p.cross_family],
        samples=data.sample_rows(),
        flags=list(data.flags),
    )


def companion_families(store: Any, family_id: str) -> List[str]:
    """Families listed in exactly the same cohort pedigrees: one edit can hold both."""
    files = family_cohort_files(store, family_id)
    out: Set[str] = set()
    for cohort in store.cohorts.values():
        if cohort.pedigree_file in files:
            out.update(cohort.families)
    out.discard(family_id)
    return sorted(f for f in out if family_cohort_files(store, f) == files)


async def open_pedigree_editor(
    anchor: ui.element,
    store: Any,
    family_id: str,
    *,
    flags: Optional[Dict[str, List[FlagRow]]] = None,
    qc: Optional[QcContext] = None,
    on_saved: Optional[Callable[[EditRecord], Any]] = None,
    on_merged: Optional[Callable[[str], Any]] = None,
) -> None:
    """Read the family's pedigree off the loop, then open the editor in ``anchor``.

    ``flags`` are the open Sample QC flags by member, for the trees' outlines.
    ``qc`` holds the family's pairs, samples, flags and cross-family pairs: the
    draft's fit and the suggestions. ``on_saved`` gets the (first family's)
    edit record once the files are written; ``on_merged`` the merged family's
    ID after a merge (without it, the editor offers none).
    """
    files = family_cohort_files(store, family_id)
    try:
        source, history = await asyncio.to_thread(
            _load, store.data_dir, family_id, files
        )
        companions = await asyncio.to_thread(companion_families, store, family_id)
    except PedigreeEditError as e:
        ui.notify(str(e), type="negative", multi_line=True)
        return
    except Exception:
        logger.exception("Failed to read the pedigree of %s", family_id)
        ui.notify(f"The pedigree of {family_id} could not be read", type="negative")
        return

    async def load_other(fid: str) -> Optional[Tuple[FamilySource, QcContext]]:
        try:
            return await asyncio.to_thread(_other_family, store, fid)
        except PedigreeEditError as e:
            ui.notify(str(e), type="negative", multi_line=True)
        except Exception:
            logger.exception("Failed to read the pedigree of %s", fid)
            ui.notify(f"The pedigree of {fid} could not be read", type="negative")
        return None

    qc = qc or QcContext()

    async def start_merge(others: List[str]) -> None:
        await open_merge_editor(
            anchor, store, [family_id, *others], on_merged=on_merged
        )

    anchor.clear()  # an earlier editor's dialog
    with anchor:
        _editor_dialog(
            Side(source, source.draft(), qc, flags or {}),
            history,
            companions,
            load_other,
            on_saved,
            start_merge=start_merge if on_merged is not None else None,
        ).open()


def _editor_dialog(
    main: Side,
    history: List[EditRecord],
    companions: List[str],
    load_other: Callable[[str], Any],
    on_saved: Optional[Callable[[EditRecord], Any]],
    *,
    start_merge: Optional[Callable[[List[str]], Any]] = None,
    merge: Optional[MergeSource] = None,
    curation: Optional[Dict[Path, int]] = None,
    on_merged: Optional[Callable[[str], Any]] = None,
) -> ui.dialog:
    family_id = main.family_id
    sides: List[Side] = [main]
    state: Dict[str, Any] = {
        "undo": [],
        "answered": set(),
        "target": merge.target if merge is not None else family_id,
    }

    def drafts() -> List[PedigreeDraft]:
        return [s.draft for s in sides]

    def snapshot() -> List[PedigreeDraft]:
        return [d.copy() for d in drafts()]

    def restore_drafts(copies: List[PedigreeDraft]) -> None:
        for side, d in zip(sides, copies):
            side.draft = d

    def change_all(change: Callable[..., Any]) -> None:
        """Run ``change(*drafts)``; undone whole on an error."""
        before = snapshot()
        try:
            change(*drafts())
        except PedigreeEditError as e:
            restore_drafts(before)
            ui.notify(str(e), type="warning")
        else:
            state["undo"].append(before)
        refresh()

    def on_side(index: int) -> Callable[[Callable[[PedigreeDraft], Any]], None]:
        def apply(change: Callable[[PedigreeDraft], Any]) -> None:
            change_all(lambda *ds: change(ds[index]))

        return apply

    def both_qc() -> Optional[QcContext]:
        if len(sides) < 2:
            return None
        return main.qc.with_family(sides[1].qc, family_id, sides[1].family_id)

    def fit_of(which: str) -> Tuple[Optional[Fit], Optional[Fit]]:
        """The fit now and with the changes, over every family open."""
        if len(sides) == 1:
            if not main.qc.pairs:
                return None, None
            return main.qc.fit(main.original), main.qc.fit(main.draft)
        qc = both_qc()
        now = qc.fit(combined_draft([s.original for s in sides]))  # type: ignore[union-attr]
        after = qc.fit(combined_draft(drafts()))  # type: ignore[union-attr]
        return now, after

    with (
        ui.dialog() as dialog,
        ui.card()
        .classes("p-4 gap-3")
        .style("width: min(1400px, 96vw); max-width: 96vw")
        .mark("pedigree-editor"),
    ):
        with ui.row().classes("w-full items-center gap-2"):
            title = ui.label(f"Edit the pedigree of family {family_id}").classes(
                "text-lg font-semibold text-blue-800"
            )
            if merge is not None:

                def on_target(e: Any) -> None:
                    state["target"] = (e.value or "").strip()
                    main.draft.family_id = state["target"]
                    update_title()
                    summary.refresh()

                ui.input(
                    "Merged family ID", value=merge.target, on_change=on_target
                ).props("dense outlined").classes("w-64").mark("ped-merge-fid")
            writes = ui.label("").classes("text-xs text-gray-500")
            ui.space()
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")

        def update_title() -> None:
            if merge is not None:
                title.set_text(f"Merge families {', '.join(merge.family_ids)}")
                writes.set_text(
                    "Writes "
                    + ", ".join(f"{c}.pedigree.tsv" for c in merge.cohorts)
                    + f" and {state['target']}.pedigree.tsv"
                )
                return
            fids = [s.family_id for s in sides]
            title.set_text(
                f"Edit the pedigree of family {fids[0]}"
                if len(fids) == 1
                else f"Edit the pedigrees of families {' and '.join(fids)}"
            )
            writes.set_text(
                "Writes "
                + ", ".join(f"{c}.pedigree.tsv" for c in main.source.cohorts)
                + " and "
                + ", ".join(f"{f}.pedigree.tsv" for f in fids)
            )

        update_title()

        @ui.refreshable
        def trees() -> None:
            for index, side in enumerate(sides):
                changed = changed_ids(side.original, side.draft)
                tag = "" if index == 0 else "2"
                if len(sides) > 1:
                    ui.label(f"Family {side.family_id}").classes(
                        "text-sm font-semibold text-blue-900"
                    )
                with ui.row().classes("w-full no-wrap gap-4 items-start"):
                    for heading, which, kind in (
                        ("Now", side.original, "now"),
                        ("After your changes", side.draft, "new"),
                    ):
                        with ui.column().classes(
                            "flex-1 min-w-0 border rounded p-2 bg-gray-50 gap-1"
                        ):
                            ui.label(heading).classes(
                                "text-sm font-semibold text-gray-700"
                            )
                            ui.html(
                                draft_tree_html(
                                    which,
                                    f"ped-edit{tag}-{kind}-",
                                    side.flags,
                                    changed if kind == "new" else (),
                                ),
                                sanitize=False,
                            ).classes("w-full").mark(f"ped-edit{tag}-{kind}")
            ui.label(
                "Orange outline: changed or added. Dashed: not sequenced."
                " Coloured outline: open Sample QC issue."
            ).classes("text-xs text-gray-500")

        trees()

        @ui.refreshable
        def suggestions() -> None:
            found = (
                [(0, s) for s in main.qc.suggest(main.draft)]
                if (main.qc.pairs or main.qc.samples)
                else []
            )
            if len(sides) > 1:
                other = sides[1]
                if other.qc.pairs or other.qc.samples:
                    found += [(1, s) for s in other.qc.suggest(other.draft)]
            across = (
                suggest_across(main.draft, sides[1].draft, both_qc())  # type: ignore[arg-type]
                if len(sides) > 1
                else []
            )
            if not found and not across and not (main.qc.pairs or main.qc.samples):
                return
            with (
                ui.column()
                .classes("w-full gap-1 border rounded p-2 bg-blue-50")
                .mark("ped-suggestions")
            ):
                ui.label("Suggestions from Sample QC").classes(
                    "text-sm font-semibold text-blue-900"
                )
                if not found and not across:
                    ui.label(
                        "None: Sample QC suggests no other change for this pedigree."
                    ).classes("text-sm text-gray-600 italic")
                for c in across:
                    _suggestion_row(c, lambda c=c: take_across(c))
                for index, s in found:
                    prefix = f"{sides[index].family_id}: " if len(sides) > 1 else ""
                    _suggestion_row(s, lambda s=s, i=index: take(s, i), prefix)

        suggestions()

        @ui.refreshable
        def companion_picker() -> None:
            if merge is not None or len(sides) > 1 or not companions:
                return
            related = [f for f in main.qc.related_families() if f in companions]
            options = {f: f"{f} (related in Sample QC)" for f in related}
            options.update({f: f for f in companions if f not in options})
            with ui.row().classes("w-full items-center gap-2"):
                ui.label("Another family:").classes("text-sm text-gray-700")
                pick = (
                    ui.select(
                        options,
                        value=related[0] if related else None,
                        with_input=True,
                    )
                    .props("dense outlined options-dense")
                    .classes("min-w-[14rem]")
                    .mark("ped-other-family")
                )

                async def open_other() -> None:
                    fid = pick.value
                    if not fid:
                        return
                    loaded = await load_other(fid)
                    if loaded is None:
                        return
                    source, qc = loaded
                    sides.append(
                        Side(source, source.draft(), qc, flags_by_sample(qc.flags))
                    )
                    state["undo"].clear()
                    update_title()
                    companion_picker.refresh()
                    refresh()

                ui.button(
                    "Open alongside", icon="call_split", on_click=open_other
                ).props("outline dense").mark("ped-other-open")
                ui.label(
                    "To move a sample to it or swap one with it: both families"
                    " then need ghfc-ngs --clean-stale-families."
                ).classes("text-xs text-gray-500")
            if start_merge is None:
                return
            with ui.row().classes("w-full items-center gap-2"):
                ui.label("Merge with:").classes("text-sm text-gray-700")
                merging = (
                    ui.select(
                        options,
                        value=related[:1],
                        multiple=True,
                        with_input=True,
                    )
                    .props("dense outlined options-dense use-chips")
                    .classes("min-w-[14rem]")
                    .mark("ped-merge-with")
                )

                async def begin_merge() -> None:
                    chosen = list(merging.value or [])
                    if not chosen:
                        ui.notify("Choose the families to merge with", type="warning")
                        return
                    dialog.close()
                    await start_merge(chosen)

                ui.button("Merge…", icon="merge", on_click=begin_merge).props(
                    "outline dense"
                ).mark("ped-merge-start")
                ui.label(
                    "One family under a new ID (or one of theirs); their"
                    " curation follows."
                ).classes("text-xs text-gray-500")

        companion_picker()

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
            for index, side in enumerate(sides):
                d = side.draft
                changed = changed_ids(side.original, d)
                other = sides[1 - index] if len(sides) > 1 else None
                if len(sides) > 1:
                    ui.label(f"Family {side.family_id}").classes(
                        "text-sm font-semibold text-blue-900 pt-2"
                    )
                with ui.grid(
                    columns="minmax(9rem, 1.4fr) repeat(4, minmax(8rem, 1fr)) auto"
                ).classes("w-full items-center gap-x-3 gap-y-1"):
                    for header in (
                        "Member",
                        "Father",
                        "Mother",
                        "Sex",
                        "Phenotype",
                        "",
                    ):
                        ui.label(header).classes("text-xs font-semibold text-gray-500")
                    for m in [*d.members, *d.implied]:
                        _member_row(
                            d,
                            m,
                            m.iid in changed,
                            m not in d.members,
                            on_side(index),
                            other.draft if other is not None else None,
                            lambda change, i=index: change_all(
                                (lambda a, b: change(a, b))
                                if i == 0
                                else (lambda a, b: change(b, a))
                            ),
                        )
                _shared_dummy_controls(
                    d, on_side(index), "ped" if index == 0 else "ped2"
                )

        rows()

        @ui.refreshable
        def summary() -> None:
            originals = [s.original for s in sides]
            if merge is not None:
                target = state["target"]
                changes = {main.family_id: merge_changes(merge, main.draft, target)}
                any_change = True
                found = merge_problems(merge, main.draft, target)
            else:
                changes = {
                    s.family_id: describe_changes(s.original, s.draft) for s in sides
                }
                any_change = any(changes.values())
                found = cross_problems(
                    originals, drafts(), outside_ids([s.source for s in sides])
                )
            with ui.row().classes("w-full no-wrap gap-4 items-start"):
                with ui.column().classes("flex-1 min-w-0 gap-1").mark("ped-changes"):
                    ui.label("Changes").classes("text-sm font-semibold text-gray-700")
                    if not any_change:
                        ui.label("None yet.").classes("text-sm text-gray-500 italic")
                    for side in sides:
                        clauses = changes[side.family_id]
                        if not clauses:
                            continue
                        if len(sides) > 1:
                            ui.label(f"Family {side.family_id}").classes(
                                "text-xs font-semibold text-gray-600"
                            )
                        for clause in clauses:
                            ui.label(clause).classes("text-sm font-mono")
                        if merge is not None:
                            text = merge_rerun_text(merge, state["target"])
                        else:
                            text = RERUN_TEXT[rerun_needed(side.original, side.draft)]
                        ui.label(text).classes("text-xs text-blue-800").mark(
                            f"ped-rerun-{side.family_id}"
                        )
                    if merge is not None and curation:
                        ui.label(
                            "Curation moved to the merged family: "
                            + ", ".join(
                                f"{n} in {p.relative_to(merge.data_dir)}"
                                for p, n in curation.items()
                            )
                            + ". Each file is copied into its history first."
                        ).classes("text-xs text-gray-700").mark("ped-merge-curation")
                    now, after = fit_of("all")
                    if now is not None:
                        ui.label(fit_text(now, after, any_change)).classes(
                            "text-xs text-gray-700"
                        ).mark("ped-fit")
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
            state["blocked"] = not any_change or any(p.level == ERROR for p in found)
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

        if history and merge is None:
            _history(main.source, history, main.flags, lambda text: restore(text))

    def update_save() -> None:
        ready = not state.get("blocked", True) and bool((comment.value or "").strip())
        save_button.set_enabled(ready)

    state["update_save"] = update_save

    def refresh() -> None:
        trees.refresh()
        suggestions.refresh()
        rows.refresh()
        summary.refresh()

    def take(s: Suggestion, index: int = 0) -> None:
        state["answered"].update(s.flags)
        on_side(index)(s.apply)

    def take_across(c: CrossSuggestion) -> None:
        state["answered"].update(c.flags)
        change_all(lambda a, b: c.apply(a, b))

    def undo() -> None:
        if state["undo"]:
            restore_drafts(state["undo"].pop())
            refresh()

    def reset() -> None:
        if any(s.draft.members != s.original.members for s in sides):
            state["undo"].append(snapshot())
            for s in sides:
                s.draft = s.source.draft()
            state["answered"].clear()
            refresh()

    def restore(text: str) -> None:
        def change(d: PedigreeDraft) -> None:
            restored = draft_from_text(text, family_id)
            d.members = restored.members

        on_side(0)(change)
        ui.notify("Rows restored in the draft: review them, then save with a comment")

    async def save() -> None:
        save_button.disable()
        if merge is not None:
            await save_merged()
            return
        try:
            if len(sides) == 1:
                records = [
                    await save_pedigree_edit_async(
                        main.source,
                        main.draft,
                        user=get_current_user(),
                        comment=comment.value or "",
                        flags=sorted(state["answered"]),
                    )
                ]
            else:
                records = await save_cross_family_edit_async(
                    [(s.source, s.draft) for s in sides],
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
        reruns = {r.rerun for r in records}
        rerun = RERUN_CLEAN if RERUN_CLEAN in reruns else records[0].rerun
        ui.notify(
            f"Pedigree of {' and '.join(r.family_id for r in records)} saved. "
            + RERUN_TEXT[rerun],
            type="positive",
            multi_line=True,
        )
        dialog.close()
        if on_saved is not None:
            result = on_saved(records[0])
            if asyncio.iscoroutine(result):
                await result

    async def save_merged() -> None:
        target = state["target"]
        try:
            await save_merge_async(
                merge,  # type: ignore[arg-type]
                main.draft,
                target,
                user=get_current_user(),
                comment=comment.value or "",
                flags=sorted(state["answered"]),
            )
        except (PermissionError, PedigreeEditError) as e:
            ui.notify(str(e), type="negative", multi_line=True)
            update_save()
            return
        except Exception:
            logger.exception("Failed to merge %s", ", ".join(merge.family_ids))  # type: ignore[union-attr]
            ui.notify("The families could not be merged", type="negative")
            update_save()
            return
        ui.notify(
            f"Families merged into {target}. " + merge_rerun_text(merge, target),  # type: ignore[arg-type]
            type="positive",
            multi_line=True,
        )
        dialog.close()
        if on_merged is not None:
            result = on_merged(target)
            if asyncio.iscoroutine(result):
                await result

    update_save()
    return dialog


def merge_rerun_text(merge: MergeSource, target: str) -> str:
    """What the next ghfc-ngs run does for a merge, worded for the curator."""
    if merge_rerun(merge, target) == RERUN_CLEAN:
        return (
            f"{target} gains members: ghfc-ngs needs --clean-stale-families. The"
            " other families' directories stay, orphaned (listed on the Workflow"
            " tab)."
        )
    return (
        f"{target} is a new family for ghfc-ngs: a normal run builds it from the"
        " members' gVCFs. The old families' directories stay, orphaned (listed on"
        " the Workflow tab)."
    )


def _load_merge(
    store: Any, family_ids: List[str]
) -> Tuple[MergeSource, QcContext, Dict[Path, int], Dict[str, List[FlagRow]]]:
    """Blocking: the families to merge, their QC as one, the curation to
    re-key, and the open flags by member (for the trees' outlines)."""
    files = [family_cohort_files(store, fid) for fid in family_ids]
    merge = load_merge(store.data_dir, family_ids, files)
    contexts: Dict[str, QcContext] = {}
    label = None
    for fid in family_ids:
        data = load_family_qc(store.data_dir, fid)
        if data is None:
            contexts[fid] = QcContext()
            continue
        label = label or data.label
        contexts[fid] = QcContext(
            pairs=[p for p in pair_rows(data.pairs) if not p.cross_family],
            samples=data.sample_rows(),
            flags=list(data.flags),
        )
    cohort_pairs: List[PairRow] = []
    for cohort in store.cohorts.values():
        if cohort.pedigree_file in files[0]:
            qc = load_cohort_qc(cohort, label)
            if qc is not None:
                cohort_pairs += pair_rows(qc.pairs)
    qc = merged_context(contexts, cohort_pairs, merge.target)
    acks = load_acks(store.data_dir)
    open_flags = flags_by_sample(f for f in qc.flags if f not in acks)
    return merge, qc, curation_counts(store.data_dir, family_ids), open_flags


async def open_merge_editor(
    anchor: ui.element,
    store: Any,
    family_ids: List[str],
    *,
    on_merged: Optional[Callable[[str], Any]] = None,
) -> None:
    """Read the families off the loop, then open the editor on them merged.

    The merged family is one draft under the suggested ID (the IDs joined by
    ``_``, editable), with the families' Sample QC read as one, so its members
    can be linked before saving. ``on_merged`` gets the merged family's ID.
    """
    try:
        merge, qc, curation, open_flags = await asyncio.to_thread(
            _load_merge, store, family_ids
        )
    except PedigreeEditError as e:
        ui.notify(str(e), type="negative", multi_line=True)
        return
    except Exception:
        logger.exception("Failed to read the families %s", family_ids)
        ui.notify("The families could not be read", type="negative")
        return
    source = merge.as_source()
    anchor.clear()
    with anchor:
        _editor_dialog(
            Side(source, source.draft(), qc, open_flags),
            [],
            [],
            lambda fid: None,
            None,
            merge=merge,
            curation=curation,
            on_merged=on_merged,
        ).open()


def render_merge_candidates(
    anchor: ui.element,
    store: Any,
    cross_pairs: List[PairRow],
    on_merged: Callable[[str], Any],
    limit: int = 20,
) -> None:
    """Curators, on the cohort's Sample QC tab: groups of families joined by
    first- or second-degree relatives across them, each with Merge….

    A group whose families are not in the same cohort pedigrees cannot be
    merged in one edit: its button is disabled, with the reason.
    """
    groups = merge_candidates(cross_pairs)
    if not groups:
        ui.label(
            "None: no first- or second-degree relatives across the families shown."
        ).classes("text-xs text-gray-500")
        return
    ui.label(
        f"{len(groups)} group(s) of families whose members are close relatives:"
        " one family entered under several IDs, or families that should be one."
    ).classes("text-xs text-gray-600")
    with ui.column().classes("w-full gap-1").mark("merge-candidates"):
        for group in groups[:limit]:
            fids = list(group.families)
            files = [family_cohort_files(store, f) for f in fids]
            together = all(f == files[0] for f in files[1:])
            with ui.row().classes("w-full items-start no-wrap gap-2 border-b py-1"):
                with ui.column().classes("flex-1 min-w-0 gap-0"):
                    ui.label(" + ".join(fids)).classes("text-sm font-semibold")
                    for line in group.evidence()[:3]:
                        ui.label(line).classes("text-xs text-gray-600")
                    if len(group.pairs) > 3:
                        ui.label(f"and {len(group.pairs) - 3} more pair(s)").classes(
                            "text-xs text-gray-500"
                        )
                button = (
                    ui.button(
                        "Merge…",
                        icon="merge",
                        on_click=lambda f=fids: open_merge_editor(
                            anchor, store, f, on_merged=on_merged
                        ),
                    )
                    .props("outline dense size=sm")
                    .mark(f"merge-candidate-{'+'.join(fids)}")
                )
                if not together:
                    button.disable()
                    button.tooltip("Not in the same cohort pedigrees: merge by hand")
        if len(groups) > limit:
            ui.label(f"… and {len(groups) - limit} more group(s)").classes(
                "text-xs text-gray-500"
            )


def _member_row(
    d: PedigreeDraft,
    m: Member,
    changed: bool,
    implied: bool,
    apply: Callable[[Callable[[PedigreeDraft], Any]], None],
    other: Optional[PedigreeDraft] = None,
    across: Optional[
        Callable[[Callable[[PedigreeDraft, PedigreeDraft], Any]], None]
    ] = None,
) -> None:
    """One member's row. ``other`` is the family open alongside; ``across``
    runs a change on (this family's draft, the other's)."""
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
                .tooltip("Swap or move: for a sample swap")
                .mark(f"ped-swap-{iid}")
            ):
                with ui.menu():
                    for o in others:
                        ui.menu_item(
                            f"Swap {iid} and {o}",
                            on_click=lambda o=o: apply(lambda dr: dr.swap(iid, o)),
                        ).mark(f"ped-swap-{iid}-{o}")
                    if other is not None and across is not None:
                        ui.separator()
                        for o in other.sequenced_ids:
                            ui.menu_item(
                                f"Swap {iid} and {o} ({other.family_id})",
                                on_click=lambda o=o: across(
                                    lambda mine, theirs, o=o: swap_across(
                                        mine, theirs, iid, o
                                    )
                                ),
                            ).mark(f"ped-swap-{iid}-{o}")
                        ui.menu_item(
                            f"Move {iid} to {other.family_id}",
                            on_click=lambda: across(
                                lambda mine, theirs: move_member(mine, theirs, iid)
                            ),
                        ).mark(f"ped-move-{iid}")
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
    d: PedigreeDraft,
    apply: Callable[[Callable[[PedigreeDraft], Any]], None],
    mark: str = "ped",
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
            .mark(f"{mark}-shared-children")
        )
        role = (
            ui.select(
                {"both": "Father and mother", FATHER: "Father", MOTHER: "Mother"},
                value="both",
            )
            .props("dense outlined options-dense")
            .mark(f"{mark}-shared-role")
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
        ).mark(f"{mark}-shared-add")
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


def _suggestion_row(
    s: "Suggestion | CrossSuggestion", take: Callable[[], None], prefix: str = ""
) -> None:
    with ui.row().classes("w-full items-start no-wrap gap-2"):
        with ui.column().classes("flex-1 min-w-0 gap-0"):
            with ui.row().classes("items-center gap-2"):
                ui.label(prefix + s.title).classes("text-sm font-semibold")
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


def _waiting(cohort: Any) -> Tuple[List[WaitingFamily], List[Orphan]]:
    """Blocking: the cohort's families edited since its last run started, and
    the directories of the families merged away from it."""
    state = read_cohort_state(cohort.pedigree_file.parent)
    last = state.last_run if state is not None else None
    members = {
        fid: [s.sample_id for s in family.samples]
        for fid, family in cohort.families.items()
    }
    data_dir = cohort.pedigree_file.parents[2]
    waiting = waiting_families(
        data_dir, members, last.started_at if last is not None else None
    )
    records = [r for r in read_edit_index(data_dir) if cohort.name in r.cohorts]
    return waiting, orphaned_families(data_dir, records, cohort.families)


def _drift_text(w: WaitingFamily) -> str:
    if w.missing is None:
        return "no joint call yet"
    parts = []
    if w.missing:
        parts.append(f"the joint call lacks {', '.join(w.missing)}")
    if w.extra:
        parts.append(f"has {', '.join(w.extra)} no longer in the family")
    return "; ".join(parts) if parts else "the joint call has the same members"


def render_waiting_families(cohort: Any) -> None:
    """Curators, on the Workflow tab: families whose pedigree was edited since
    the last ghfc-ngs run started, and which need ``--clean-stale-families``.

    Read off the event loop (the edit index, then each family's joint-call
    header); nothing is drawn when no family waits.
    """
    box = ui.column().classes("w-full min-w-0 gap-1")

    async def load() -> None:
        try:
            waiting, orphans = await asyncio.to_thread(_waiting, cohort)
        except Exception:
            logger.warning(
                "Failed to list the pedigree edits of %s", cohort.name, exc_info=True
            )
            return
        if orphans:
            with (
                box,
                ui.column()
                .classes(
                    "w-full gap-1 border border-gray-300 rounded-md p-3 bg-gray-50"
                )
                .mark("workflow-orphaned-families"),
            ):
                ui.label(
                    f"Directories of merged families ({len(orphans)}): ghfc-ngs"
                    " leaves them. Remove them once the merged family is built."
                ).classes("text-sm font-semibold text-gray-800")
                for o in orphans:
                    ui.label(
                        f"{o.family_id} → {o.into} (merged {o.timestamp[:10]}): {o.path}"
                    ).classes("text-xs text-gray-700 font-mono")
        if not waiting:
            return
        clean = [w.family_id for w in waiting if w.needs_clean]
        with (
            box,
            ui.column()
            .classes("w-full gap-1 border border-amber-300 rounded-md p-3 bg-amber-50")
            .mark("workflow-pedigree-waiting"),
        ):
            ui.label(
                f"Pedigree edits waiting for a ghfc-ngs run"
                f" ({len(waiting)} famil{'y' if len(waiting) == 1 else 'ies'})"
            ).classes("text-sm font-semibold text-amber-900")
            for w in waiting:
                need = (
                    "needs --clean-stale-families" if w.needs_clean else "a normal run"
                )
                ui.label(
                    f"{w.family_id}: edited {w.last.timestamp[:10]} by"
                    f" {w.last.user} ({w.last.comment}); {need}; {_drift_text(w)}"
                ).classes("text-xs text-gray-800")
            text = "Run ghfc-ngs as usual"
            if clean:
                text += (
                    ", with --clean-stale-families for "
                    + ", ".join(clean)
                    + " (--dry-run first lists what it deletes)"
                )
            ui.label(text + ".").classes("text-xs text-amber-900 font-semibold")

    ui.timer(0.05, load, once=True)
