"""A cohort's pedigree edit history, as one self-contained HTML document.

Every edit genetics-viz saved to a family of the cohort (``Cohorts`` of
``validations/pedigree_edits.tsv``), grouped by family: the tree before the
family's first edit and after its last (the changed members outlined), then
each edit with its date, curator, comment, changes and the rerun it needs.

The trees come from the rows each edit kept (``pedigree_edits/<FID>.<EditId>
.{before,after}.tsv``), so the document says what the files held, not what
the cohort pedigree holds now. Built off the event loop (:func:`build_history`
reads files); downloaded from the cohort's Sample QC tab, printed to PDF from
the browser like the family report.
"""

import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from markupsafe import Markup

from .pedigree_edit import (
    FATHER,
    MOTHER,
    RERUN_CLEAN,
    RERUN_NONE,
    RERUN_NORMAL,
    UNKNOWN,
    EditRecord,
    PedigreeDraft,
    draft_from_text,
    edit_rows,
    read_edit_index,
    read_family_edits,
)
from .pedigree_labels import normalize_pheno, normalize_sex
from .pedigree_svg import PedigreePerson, pedigree_svg
from .report_html import render_template

logger = logging.getLogger(__name__)

#: The rerun an edit needs, worded for the document.
RERUN_WORDS = {
    RERUN_NONE: "nothing to rerun (phenotype only)",
    RERUN_NORMAL: "a normal ghfc-ngs run",
    RERUN_CLEAN: "a run with --clean-stale-families",
    "merged": "none for this ID (merged or renamed)",
}

#: A tree's box in the document: two side by side on an A4 page.
MM_PER_UNIT = 14
MAX_WIDTH_MM = 88
MAX_HEIGHT_MM = 80


def changed_ids(original: PedigreeDraft, draft: PedigreeDraft) -> Set[str]:
    """Members whose row the edit adds or changes (codes normalised)."""
    was = {m.iid: m.key() for m in [*original.members, *original.implied]}
    return {
        m.iid for m in [*draft.members, *draft.implied] if was.get(m.iid) != m.key()
    }


def _parent(code: str) -> str:
    return "" if code == UNKNOWN else code


def draft_people(
    draft: PedigreeDraft, classes: Optional[Mapping[str, Sequence[str]]] = None
) -> List[PedigreePerson]:
    """The draft's members (and implied dummies) as drawing people.

    ``classes``: extra classes per member, on its ``<g class="node …">``.
    """
    return [
        PedigreePerson(
            m.iid,
            _parent(m.parent(FATHER)),
            _parent(m.parent(MOTHER)),
            normalize_sex(m.sex),
            normalize_pheno(m.phenotype),
            dummy=m.is_dummy,
            classes=tuple((classes or {}).get(m.iid, ())),
        )
        for m in [*draft.members, *draft.implied]
    ]


@dataclass
class Tree:
    """One drawing of the document, or the reason there is none."""

    title: str
    svg: str = ""
    #: CSS width of the box (``"62mm"``).
    width: str = ""
    missing: str = ""

    @property
    def markup(self) -> Markup:
        # pypedigree builds the SVG from escaped text (pedigree_svg).
        return Markup(self.svg)


@dataclass
class FamilyHistory:
    family_id: str
    #: Oldest first.
    edits: List[EditRecord]
    before: Tree
    after: Tree
    #: Last edit not followed by a ghfc-ngs run of the cohort (``None``: unknown).
    awaiting_rerun: Optional[bool] = None

    @property
    def users(self) -> List[str]:
        return sorted({e.user for e in self.edits if e.user})


@dataclass
class CohortHistory:
    cohort_name: str
    generated_at: str
    generated_by: str
    families: List[FamilyHistory] = field(default_factory=list)

    @property
    def n_edits(self) -> int:
        return sum(len(f.edits) for f in self.families)

    @property
    def users(self) -> List[str]:
        return sorted({u for f in self.families for u in f.users})

    @property
    def document_title(self) -> str:
        return f"{self.cohort_name}_pedigree_history_{self.generated_at[:10]}"


def cohort_edits(data_dir: Path, cohort_name: str) -> List[EditRecord]:
    """The data directory's edits that wrote this cohort's pedigree, oldest first.

    Listed by the index; each record is taken from its family's own log when
    there, which also keeps the Sample QC flags the edit addressed.
    """
    listed = [r for r in read_edit_index(data_dir) if cohort_name in r.cohorts]
    logs: Dict[str, Dict[str, EditRecord]] = {}
    out = []
    for record in listed:
        if record.family_id not in logs:
            logs[record.family_id] = {
                r.edit_id: r for r in read_family_edits(data_dir, record.family_id)
            }
        out.append(logs[record.family_id].get(record.edit_id, record))
    return out


def _rows(data_dir: Path, record: EditRecord, which: str) -> Optional[PedigreeDraft]:
    try:
        text = edit_rows(data_dir, record.family_id, record.edit_id, which)
    except OSError:
        return None
    return draft_from_text(text, record.family_id)


def _tree(
    title: str,
    draft: Optional[PedigreeDraft],
    id_prefix: str,
    changed: Set[str] = frozenset(),
    empty: str = "",
) -> Tree:
    if draft is None:
        return Tree(title, missing="The rows kept for this edit are missing.")
    if not draft.members and not draft.implied:
        return Tree(title, missing=empty or "No rows under this family ID.")
    people = draft_people(draft, {iid: ("changed",) for iid in changed})
    drawing = pedigree_svg(draft.family_id, people, id_prefix=id_prefix)
    if drawing is None or not drawing.width or not drawing.height:
        return Tree(title, missing="This pedigree cannot be drawn.")
    width = min(
        drawing.width * MM_PER_UNIT,
        MAX_WIDTH_MM,
        MAX_HEIGHT_MM * drawing.width / drawing.height,
    )
    return Tree(title, svg=drawing.svg, width=f"{width:.0f}mm")


def _merged_note(record: EditRecord) -> str:
    """What became of a family merged into another, from its record."""
    return next(iter(record.changes), "") or "Merged into another family."


def _family(
    data_dir: Path, family_id: str, edits: List[EditRecord], index: int
) -> FamilyHistory:
    first, last = edits[0], edits[-1]
    before = _rows(data_dir, first, "before")
    after = _rows(data_dir, last, "after")
    changed = changed_ids(before, after) if before and after else set()
    prefix = f"h{index}-"
    return FamilyHistory(
        family_id=family_id,
        edits=edits,
        before=_tree(
            "Before the first edit",
            before,
            prefix + "b-",
            empty="No rows under this ID before: the family was created by a merge"
            " or a rename.",
        ),
        after=_tree(
            "After the last edit" if len(edits) > 1 else "After the edit",
            after,
            prefix + "a-",
            changed,
            empty=_merged_note(last),
        ),
    )


def build_history(
    data_dir: Path,
    cohort_name: str,
    *,
    user: str = "",
    statuses: Optional[Mapping[str, Any]] = None,
    now: Optional[datetime] = None,
) -> CohortHistory:
    """Blocking: the cohort's edits grouped by family, each with its trees.

    ``statuses``: ``FID -> EditStatus`` (``family_status.cohort_edits``), for
    whether each family's last edit still waits for a run.
    """
    by_family: Dict[str, List[EditRecord]] = {}
    for record in cohort_edits(data_dir, cohort_name):
        by_family.setdefault(record.family_id, []).append(record)
    # Most recently edited family first: the document is read for what changed.
    order: List[Tuple[str, List[EditRecord]]] = sorted(
        by_family.items(), key=lambda kv: (kv[1][-1].timestamp, kv[0]), reverse=True
    )
    families = []
    for index, (fid, edits) in enumerate(order):
        family = _family(data_dir, fid, edits, index)
        status = (statuses or {}).get(fid)
        if status is not None:
            family.awaiting_rerun = bool(status.awaiting_rerun)
        families.append(family)
    moment = now or datetime.now()
    return CohortHistory(
        cohort_name=cohort_name,
        generated_at=moment.strftime("%Y-%m-%d %H:%M"),
        generated_by=user,
        families=families,
    )


def render_history_html(history: CohortHistory) -> str:
    """The whole document, ``<!DOCTYPE html>`` included."""
    return render_template(
        "pedigree_history.html",
        history=history,
        rerun_words=RERUN_WORDS,
        when=lambda ts: ts.replace("T", " "),
    )


def history_html(
    data_dir: Path,
    cohort_name: str,
    *,
    user: str = "",
    statuses: Optional[Mapping[str, Any]] = None,
) -> str:
    """Blocking: the document, or ``""`` when no edit wrote this cohort."""
    history = build_history(data_dir, cohort_name, user=user, statuses=statuses)
    return render_history_html(history) if history.families else ""
