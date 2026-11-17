"""Curator edits of a family's pedigree, written so ghfc-ngs reruns the family.

ghfc-ngs reads one pedigree, the cohort's ``<NAME>.pedigree.tsv``, and writes
each family's ``<FID>.pedigree.tsv`` from it (``FAMILIAL_PEDIGREE``: the
cohort's header line, then the family's rows verbatim). An edit is therefore
made in the cohort file, in place: every line of another family stays byte for
byte, and the family file is written as the pipeline would write it, so the
next run finds it up to date. Every step that reads the pedigree keeps an
``inputs.tsv`` snapshot of the family's rows (``lib/InputsSnapshot.groovy``)
and reruns when they change, so a corrected parent or sex needs a normal run;
only a change of sequenced members needs ``--clean-stale-families``, and this
module refuses one for now.

The model is a :class:`PedigreeDraft`: the family's rows, edited by a handful
of operations (parents, sex, phenotype, dummy parents, a swap). :func:`problems`
checks it, :func:`describe_changes` words it as the pipeline's run summary does,
and :func:`save_pedigree_edit` writes it. Each save keeps the family's rows
before and after under ``families/…/<FID>/pedigree_edits/`` with a log line
(who, when, why, what), and one line in ``validations/pedigree_edits.tsv`` for
the cohort page. The plan is the vault note ``Pedigree editor``.
"""

from __future__ import annotations

import asyncio
import csv
import fcntl
import hashlib
import os
import re
import shutil
import tempfile
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from genetics_viz.models import DUMMY_PREFIX, is_dummy_id
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.sharding import get_family_path
from genetics_viz.utils.tsv import append_tsv_row

#: The codes the pipeline reads as an unknown parent (``InputsSnapshot.pedigreeRow``).
MISSING_PARENT_CODES = frozenset({"", "0", "-9", "."})
#: What the editor writes for an unknown parent, sex or phenotype.
UNKNOWN = "0"
MALE, FEMALE = "1", "2"
UNAFFECTED, AFFECTED = "1", "2"

FATHER, MOTHER = "father", "mother"
#: The letter a dummy parent's ID ends with, per role (``lib/DummyMember.groovy``).
ROLE_LETTER = {FATHER: "F", MOTHER: "M"}
ROLE_SEX = {FATHER: MALE, MOTHER: FEMALE}

#: A dummy added by the editor: ``DUMMY_<child>_<path>``.
NEW_DUMMY_ID = re.compile(rf"^{DUMMY_PREFIX}\S+_[FM]+$")

#: Present while a ghfc-ngs run of the cohort is live (vault ``Running the pipeline``).
PROGRESS_FILENAME = ".ghfc-ngs.progress.json"
EDITS_DIRNAME = "pedigree_edits"
LOCK_FILENAME = ".pedigree_edit.lock"

LOG_HEADER = [
    "EditId", "Timestamp", "User", "Comment", "Changes", "Cohorts", "Rerun", "Flags",
]  # fmt: skip
INDEX_HEADER = [
    "FID", "EditId", "Timestamp", "User", "Comment", "Changes", "Cohorts", "Rerun",
]  # fmt: skip

#: What the next ghfc-ngs run has to do for an edit.
RERUN_NONE = "none"  # phenotype only: no step reads it
RERUN_NORMAL = "normal"  # the steps' inputs snapshots change
RERUN_CLEAN = "clean_stale_families"  # the family's sequenced members change

ERROR, WARNING = "error", "warning"

#: Missing-parent spellings the app accepts (``models.MISSING_PARENT``) but
#: ghfc-ngs does not: it reads them as a parent's ID.
_FLOAT_MISSING = frozenset({"0.0", "-9.0"})


class PedigreeEditError(Exception):
    """An edit that cannot be made or saved; the message is for the curator."""


# ---------------------------------------------------------------------------
# Codes, as ghfc-ngs normalises them
# ---------------------------------------------------------------------------


def norm_parent(code: str) -> str:
    """An unknown parent is ``0``, whichever way it is written."""
    code = code.strip()
    return UNKNOWN if code in MISSING_PARENT_CODES else code


def _norm_12(code: str) -> str:
    code = re.sub(r"\.0+$", "", code.strip())
    return code if code in ("1", "2") else UNKNOWN


def norm_sex(code: str) -> str:
    """``1`` or ``2`` (``1.0`` and ``2.0`` read as such), else ``0``."""
    return _norm_12(code)


def norm_phenotype(code: str) -> str:
    """``1`` unaffected, ``2`` affected, else ``0`` (unknown)."""
    return _norm_12(code)


def snapshot_row(cells: Sequence[str]) -> Tuple[str, str, str, str, str]:
    """``InputsSnapshot.pedigreeRow``: fid, sample, father, mother, sex, normalised."""
    return (
        cells[0].strip(),
        cells[1].strip(),
        norm_parent(cells[2]),
        norm_parent(cells[3]),
        norm_sex(cells[4]),
    )


# ---------------------------------------------------------------------------
# The draft
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Member:
    """One pedigree row, its cells kept as written until edited."""

    iid: str
    father: str = UNKNOWN
    mother: str = UNKNOWN
    sex: str = UNKNOWN
    phenotype: str = UNKNOWN
    #: Columns after the sixth, written back unchanged.
    extra: Tuple[str, ...] = ()

    @property
    def is_dummy(self) -> bool:
        return is_dummy_id(self.iid)

    def parent(self, role: str) -> str:
        return norm_parent(self.father if role == FATHER else self.mother)

    def key(self) -> Tuple[str, str, str, str]:
        """What the pipeline's steps read: parents and sex, normalised."""
        return (
            norm_parent(self.father),
            norm_parent(self.mother),
            norm_sex(self.sex),
            norm_phenotype(self.phenotype),
        )


@dataclass
class PedigreeDraft:
    """A family's pedigree rows being edited, in file order.

    A dummy named in a father or mother column with no row of its own is not a
    member here but :attr:`implied`: it is written only once edited.
    """

    family_id: str
    members: List[Member] = field(default_factory=list)

    def copy(self) -> "PedigreeDraft":
        return PedigreeDraft(self.family_id, list(self.members))

    @property
    def ids(self) -> List[str]:
        return [m.iid for m in self.members]

    @property
    def sequenced_ids(self) -> List[str]:
        return [m.iid for m in self.members if not m.is_dummy]

    def get(self, iid: str) -> Optional[Member]:
        return next((m for m in self.members if m.iid == iid), None)

    def __contains__(self, iid: str) -> bool:
        return self.get(iid) is not None

    @property
    def implied(self) -> List[Member]:
        """Dummies named as a parent with no row: founders, sex from the role."""
        out: Dict[str, Member] = {}
        for m in self.members:
            for role in (FATHER, MOTHER):
                pid = m.parent(role)
                if is_dummy_id(pid) and pid not in self and pid not in out:
                    out[pid] = Member(pid, sex=ROLE_SEX[role])
        return list(out.values())

    def children(self, iid: str) -> List[str]:
        return [
            m.iid for m in self.members if iid in (m.parent(FATHER), m.parent(MOTHER))
        ]

    # --- Operations ------------------------------------------------------

    def _member(self, iid: str) -> Member:
        member = self.get(iid)
        if member is not None:
            return member
        implied = next((m for m in self.implied if m.iid == iid), None)
        if implied is None:
            raise PedigreeEditError(f"{iid} is not a member of family {self.family_id}")
        # Edited, an implied dummy gets the row it lacked.
        self.members.append(implied)
        return implied

    def _replace(self, iid: str, **changes: str) -> None:
        member = self._member(iid)
        index = self.ids.index(iid)
        self.members[index] = replace(member, **changes)

    def set_parent(self, child: str, role: str, parent: Optional[str]) -> None:
        """Set ``child``'s father or mother; ``None`` or ``0`` makes it unknown."""
        if role not in ROLE_LETTER:
            raise ValueError(f"Unknown parent role: {role!r}")
        parent = norm_parent(parent or "")
        if parent == child:
            raise PedigreeEditError(f"{child} cannot be their own {role}")
        self._replace(child, **{role: parent})

    def set_sex(self, iid: str, sex: str) -> None:
        self._replace(iid, sex=norm_sex(sex))

    def set_phenotype(self, iid: str, phenotype: str) -> None:
        self._replace(iid, phenotype=norm_phenotype(phenotype))

    def add_dummy(
        self, iid: str, sex: str, father: str = UNKNOWN, mother: str = UNKNOWN
    ) -> None:
        """Add a dummy row: a founder unless parents are given."""
        if not is_dummy_id(iid):
            raise PedigreeEditError(f"{iid} is not a dummy ID ({DUMMY_PREFIX}…)")
        if iid in self:
            raise PedigreeEditError(f"{iid} is already in family {self.family_id}")
        self.members.append(
            Member(
                iid, norm_parent(father), norm_parent(mother), norm_sex(sex), UNKNOWN
            )
        )

    def add_dummy_parents(
        self, children: Sequence[str], roles: Iterable[str] = (FATHER, MOTHER)
    ) -> List[str]:
        """Give ``children`` a shared dummy for each unknown parent of ``roles``.

        The dummy hangs from the first child (``DUMMY_<child>_F``, or a letter
        added to a dummy child's own path). Full sibs share one; for half sibs,
        call this once per child. A role a child already has a parent for is
        left alone for that child. Returns the dummies added.
        """
        if not children:
            return []
        for child in children:
            self._member(child)
        added: List[str] = []
        for role in roles:
            missing = [c for c in children if self._member(c).parent(role) == UNKNOWN]
            if not missing:
                continue
            dummy = dummy_parent_id(missing[0], role)
            if dummy not in self:
                self.add_dummy(dummy, ROLE_SEX[role])
                added.append(dummy)
            for child in missing:
                self.set_parent(child, role, dummy)
        return added

    def swap(self, a: str, b: str) -> None:
        """Exchange two sequenced members' places: rows and parent references.

        For a sample swap: the DNA labelled ``a`` belongs where ``b`` is. Each
        barcode takes the other's parents, children, sex and phenotype; the
        data stay keyed by barcode, so nothing else moves.
        """
        for iid in (a, b):
            member = self.get(iid)
            if member is None or member.is_dummy:
                raise PedigreeEditError(
                    f"{iid} is not a sequenced member of the family"
                )
        if a == b:
            return
        swap = {a: b, b: a}

        def renamed(code: str) -> str:
            return swap.get(code.strip(), code)

        self.members = [
            replace(
                m,
                iid=swap.get(m.iid, m.iid),
                father=renamed(m.father),
                mother=renamed(m.mother),
            )
            for m in self.members
        ]

    def remove_dummy(self, iid: str) -> None:
        """Remove a dummy row nobody names as a parent."""
        member = self.get(iid)
        if member is None or not member.is_dummy:
            raise PedigreeEditError(f"{iid} is not a dummy row of the family")
        if self.children(iid):
            raise PedigreeEditError(
                f"{iid} is the parent of {', '.join(self.children(iid))}: "
                "set their parent first"
            )
        self.members = [m for m in self.members if m.iid != iid]


def dummy_parent_id(child: str, role: str) -> str:
    """``DUMMY_<child>_F`` / ``_M``; a dummy's own parent extends its path (``_MF``)."""
    letter = ROLE_LETTER[role]
    return (
        f"{child}{letter}" if is_dummy_id(child) else f"{DUMMY_PREFIX}{child}_{letter}"
    )


# ---------------------------------------------------------------------------
# Checks, changes, rerun
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Problem:
    level: str  # ERROR blocks the save, WARNING does not
    message: str
    members: Tuple[str, ...] = ()


def _ancestor_cycle(draft: PedigreeDraft) -> Optional[str]:
    parents = {m.iid: [m.parent(FATHER), m.parent(MOTHER)] for m in draft.members}
    state: Dict[str, int] = {}  # 1 visiting, 2 done

    def visit(iid: str) -> Optional[str]:
        if state.get(iid) == 1:
            return iid
        if state.get(iid) == 2 or iid not in parents:
            return None
        state[iid] = 1
        for pid in parents[iid]:
            if pid != UNKNOWN and (found := visit(pid)):
                return found
        state[iid] = 2
        return None

    for iid in parents:
        if found := visit(iid):
            return found
    return None


def problems(
    original: PedigreeDraft,
    draft: PedigreeDraft,
    other_ids: Optional[Set[str]] = None,
) -> List[Problem]:
    """What is wrong with ``draft``: errors block the save, warnings do not.

    ``other_ids`` are the IDs of the other families of the cohort files: the
    pipeline refuses one ID with two different rows anywhere in the file.
    """
    out: List[Problem] = []

    def err(message: str, *ids: str) -> None:
        out.append(Problem(ERROR, message, ids))

    def warn(message: str, *ids: str) -> None:
        out.append(Problem(WARNING, message, ids))

    members = {m.iid: m for m in draft.members}
    implied = {m.iid: m for m in draft.implied}
    before = set(original.ids)

    seen: Set[str] = set()
    for m in draft.members:
        if m.iid in seen:
            err(f"{m.iid} has two rows", m.iid)
        seen.add(m.iid)
        cells = (m.iid, m.father, m.mother, m.sex, m.phenotype, *m.extra)
        if not m.iid or any(re.search(r"\s", c.strip()) or "\t" in c for c in cells):
            err(f"{m.iid or 'A row'}: an ID or code holds a space or a tab", m.iid)
        if m.iid not in before:
            if m.is_dummy and not NEW_DUMMY_ID.match(m.iid):
                err(f"{m.iid}: a new dummy is named DUMMY_<child>_<F|M…>", m.iid)
            if other_ids and m.iid in other_ids:
                err(f"{m.iid} is already a member of another family", m.iid)

    for m in draft.members:
        father, mother = m.parent(FATHER), m.parent(MOTHER)
        if father != UNKNOWN and father == mother:
            err(f"{m.iid}: {father} is both father and mother", m.iid)
        for role, pid in ((FATHER, father), (MOTHER, mother)):
            if pid == UNKNOWN:
                continue
            if pid == m.iid:
                err(f"{m.iid} is their own {role}", m.iid)
            parent = members.get(pid) or implied.get(pid)
            if parent is None:
                warn(
                    f"{m.iid}'s {role} {pid} has no row: the pipeline reads "
                    "them as a founder",
                    m.iid,
                )
                continue
            wrong = FEMALE if role == FATHER else MALE
            if norm_sex(parent.sex) == wrong:
                sex = "female" if wrong == FEMALE else "male"
                err(f"{pid} is {m.iid}'s {role} but is recorded as {sex}", pid, m.iid)
        for role, code in ((FATHER, m.father), (MOTHER, m.mother)):
            if code.strip() in _FLOAT_MISSING:
                warn(
                    f"{m.iid}'s {role} is written {code.strip()}, which ghfc-ngs "
                    "reads as an ID, not as unknown: set it to unknown",
                    m.iid,
                )
    if cyclic := _ancestor_cycle(draft):
        err(f"{cyclic} is their own ancestor", cyclic)

    lost = sorted(set(original.sequenced_ids) - set(draft.sequenced_ids))
    gained = sorted(set(draft.sequenced_ids) - set(original.sequenced_ids))
    if lost or gained:
        what = "; ".join(
            part for part in (
                f"removed {', '.join(lost)}" if lost else "",
                f"added {', '.join(gained)}" if gained else "",
            ) if part
        )  # fmt: skip
        err(
            f"This changes the family's sequenced members ({what}). That needs "
            "ghfc-ngs --clean-stale-families and is not supported yet.",
            *lost, *gained,
        )  # fmt: skip

    for m in draft.members:
        if m.is_dummy and not draft.children(m.iid):
            warn(f"{m.iid} is no one's parent", m.iid)
        elif not m.is_dummy and norm_sex(m.sex) == UNKNOWN:
            warn(f"{m.iid}'s sex is unknown", m.iid)
    # Children who share one known parent and have no other: full or half sibs
    # is undecided until they get dummy parents (shared for full sibs).
    for role, other in ((FATHER, MOTHER), (MOTHER, FATHER)):
        groups: Dict[str, List[str]] = {}
        for m in draft.members:
            if m.parent(role) != UNKNOWN and m.parent(other) == UNKNOWN:
                groups.setdefault(m.parent(role), []).append(m.iid)
        for pid, kids in groups.items():
            if len(kids) > 1:
                warn(
                    f"{', '.join(kids)} share their {role} {pid} and have no "
                    f"{other}: full or half sibs is undecided. Add a dummy "
                    f"{other}, one shared by full sibs.",
                    *kids,
                )
    return out


def has_errors(found: Iterable[Problem]) -> bool:
    return any(p.level == ERROR for p in found)


_FIELDS = ("father", "mother", "sex", "phenotype")


def describe_changes(original: PedigreeDraft, draft: PedigreeDraft) -> List[str]:
    """One clause per member, worded as ghfc-ngs's ``FamilyPedigree.describe``.

    ``S3: father S1 -> DUMMY_S3_F``, ``DUMMY_S3_F: added``, ``X: removed``;
    codes normalised, so ``2.0`` to ``2`` is no change. An implied dummy that
    gets a row is "added" only if the row says more than the implied one.
    """
    was = {m.iid: m.key() for m in [*original.members, *original.implied]}
    now = {m.iid: m.key() for m in [*draft.members, *draft.implied]}
    clauses: List[str] = []
    for iid in sorted(set(was) | set(now)):
        if iid not in was:
            clauses.append(f"{iid}: added")
        elif iid not in now:
            clauses.append(f"{iid}: removed")
        else:
            diffs = [
                f"{name} {a} -> {b}"
                for name, a, b in zip(_FIELDS, was[iid], now[iid])
                if a != b
            ]
            if diffs:
                clauses.append(f"{iid}: {', '.join(diffs)}")
    return clauses


def rerun_needed(original: PedigreeDraft, draft: PedigreeDraft) -> str:
    """What the next ghfc-ngs run has to do: none, normal, or clean stale families."""
    if set(original.sequenced_ids) != set(draft.sequenced_ids):
        return RERUN_CLEAN
    was = {m.iid: m.key()[:3] for m in original.members}
    now = {m.iid: m.key()[:3] for m in draft.members}
    return RERUN_NORMAL if was != now else RERUN_NONE


# ---------------------------------------------------------------------------
# The files
# ---------------------------------------------------------------------------


def _split(line: str) -> List[str]:
    return line.rstrip("\r\n").split("\t")


def _ending(line: str) -> str:
    return line[len(line.rstrip("\r\n")) :]


@dataclass
class PedigreeText:
    """A cohort pedigree as read: its lines, and where the family's rows are."""

    path: Path
    text: str
    lines: List[str]  # with their line endings
    header: Optional[str]  # the first line, when it is the pipeline's header
    family_lines: List[int]  # indexes of the family's lines
    other_ids: Set[str]  # the IDs of every other family

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.text.encode("utf-8")).hexdigest()

    @property
    def newline(self) -> str:
        endings = [_ending(line) for line in self.lines if _ending(line)]
        return endings[0] if endings else "\n"

    @classmethod
    def read(cls, path: Path, family_id: str) -> "PedigreeText":
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines(keepends=True)
        header = None
        family_lines: List[int] = []
        other_ids: Set[str] = set()
        for index, line in enumerate(lines):
            if line.startswith("#") or not line.strip():
                continue
            cells = [c.strip() for c in _split(line)]
            # As parsePedigreeFile: a header only as the first line, FID first.
            if index == 0 and cells[0] == "FID":
                header = line
                continue
            if cells[0] == family_id:
                family_lines.append(index)
            elif len(cells) > 1:
                other_ids.add(cells[1])
        return cls(path, text, lines, header, family_lines, other_ids)

    def family_draft(self, family_id: str) -> PedigreeDraft:
        """The family's rows; a member's repeated rows count once, as in ghfc-ngs."""
        draft = PedigreeDraft(family_id)
        for index in self.family_lines:
            cells = _split(self.lines[index])
            if len(cells) < 6:
                raise PedigreeEditError(
                    f"{self.path.name}, line {index + 1}: {len(cells)} columns instead of 6"
                )
            cells = [c.strip() for c in cells[:6]] + cells[6:]
            if cells[1] in draft:
                continue
            draft.members.append(
                Member(
                    cells[1], cells[2], cells[3], cells[4], cells[5], tuple(cells[6:])
                )
            )
        return draft


def _member_line(family_id: str, member: Member, newline: str) -> str:
    cells = [
        family_id,
        member.iid,
        member.father,
        member.mother,
        member.sex,
        # ghfc-ngs drops trailing empty fields: an empty phenotype makes the row
        # five columns long, and the run stops.
        member.phenotype.strip() or UNKNOWN,
        *member.extra,
    ]
    return "\t".join(cells) + newline


def render_cohort_text(
    source: PedigreeText, original: PedigreeDraft, draft: PedigreeDraft
) -> str:
    """The cohort file with the family's rows edited, every other line as it was.

    An unchanged row keeps its line (and its repeats); a changed one is
    rewritten, its repeats dropped; a removed one is dropped; new rows follow
    the family's last line.
    """
    family_id = draft.family_id
    before = {m.iid: m for m in original.members}
    after = {m.iid: m for m in draft.members}
    newline = source.newline
    rows_at = set(source.family_lines)
    last = source.family_lines[-1] if source.family_lines else len(source.lines) - 1
    written: Set[str] = set()
    out: List[str] = []
    for index, line in enumerate(source.lines):
        if index in rows_at:
            iid = _split(line)[1].strip()
            member = after.get(iid)
            if member is not None and iid not in written:
                written.add(iid)
                if member == before.get(iid):
                    out.append(line if _ending(line) else line + newline)
                else:
                    out.append(_member_line(family_id, member, newline))
            elif member is not None and member == before.get(iid):
                out.append(line if _ending(line) else line + newline)  # a repeat
        else:
            out.append(line)
        if index == last:
            if out and not _ending(out[-1]):
                out[-1] += newline
            out.extend(
                _member_line(family_id, m, newline)
                for m in draft.members
                if m.iid not in written and m.iid not in before
            )
            written.update(after)
    return "".join(out)


def family_file_text(cohort_text: str, family_id: str) -> str:
    """The family pedigree as ``FAMILIAL_PEDIGREE`` writes it from the cohort's.

    The cohort's header line if it has one, then the family's lines, the first
    row of each member only.
    """
    lines = cohort_text.splitlines(keepends=True)
    out: List[str] = []
    if lines and _split(lines[0])[0].strip() == "FID":
        out.append(lines[0])
    seen: Set[str] = set()
    for line in lines[1 if out else 0 :]:
        if line.startswith("#") or not line.strip():
            continue
        cells = [c.strip() for c in _split(line)]
        if cells[0] == family_id and len(cells) > 1 and cells[1] not in seen:
            seen.add(cells[1])
            out.append(line if _ending(line) else line + "\n")
    return "".join(out)


@dataclass
class FamilySource:
    """A family's pedigree in every cohort file that lists it, as read for an edit.

    The save refuses to write a file whose content changed since this was read.
    """

    data_dir: Path
    family_id: str
    files: List[PedigreeText]
    original: PedigreeDraft

    def draft(self) -> PedigreeDraft:
        return self.original.copy()

    @property
    def cohorts(self) -> List[str]:
        return [_cohort_name(f.path) for f in self.files]

    @property
    def other_ids(self) -> Set[str]:
        return set().union(*(f.other_ids for f in self.files)) if self.files else set()


def _cohort_name(path: Path) -> str:
    return path.name.removesuffix(".pedigree.tsv")


def load_family(
    data_dir: Path, family_id: str, cohort_files: Sequence[Path]
) -> FamilySource:
    """Read the family from each cohort pedigree that lists it.

    The family must read the same in all of them (codes normalised): ghfc-ngs
    would otherwise rewrite its family file back and forth, and one edit could
    not be applied to both.
    """
    files = [PedigreeText.read(Path(p), family_id) for p in cohort_files]
    files = [f for f in files if f.family_lines]
    if not files:
        raise PedigreeEditError(f"Family {family_id} is in no cohort pedigree")
    drafts = [f.family_draft(family_id) for f in files]
    first = {m.iid: m.key() for m in drafts[0].members}
    for f, d in zip(files[1:], drafts[1:]):
        if {m.iid: m.key() for m in d.members} != first:
            raise PedigreeEditError(
                f"Family {family_id} differs between {files[0].path.name} and "
                f"{f.path.name}; make them agree by hand first"
            )
    return FamilySource(Path(data_dir), family_id, files, drafts[0])


def family_cohort_files(store, family_id: str) -> List[Path]:
    """The pedigree of every cohort of ``store`` that lists the family."""
    return [
        cohort.pedigree_file
        for cohort in store.cohorts.values()
        if family_id in cohort.families
    ]


# ---------------------------------------------------------------------------
# Saving
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class EditRecord:
    family_id: str
    edit_id: str
    timestamp: str
    user: str
    comment: str
    changes: Tuple[str, ...]
    cohorts: Tuple[str, ...]
    rerun: str
    flags: Tuple[str, ...] = ()


def edits_dir(data_dir: Path, family_id: str) -> Path:
    return get_family_path(data_dir, family_id) / EDITS_DIRNAME


def index_file(data_dir: Path) -> Path:
    return Path(data_dir) / "validations" / "pedigree_edits.tsv"


def _write_temp(target: Path, text: str) -> Path:
    """Write ``text`` beside ``target``, with its permissions, for a rename."""
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    with os.fdopen(fd, "w", encoding="utf-8", newline="") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    if target.exists():
        shutil.copymode(target, name)
    else:
        os.chmod(name, 0o664)
    return Path(name)


def _edit_id(directory: Path, family_id: str, now: datetime) -> str:
    base = now.strftime("%Y%m%dT%H%M%S")
    edit_id, n = base, 1
    while (directory / f"{family_id}.{edit_id}.before.tsv").exists():
        n += 1
        edit_id = f"{base}-{n}"
    return edit_id


def save_pedigree_edit(
    source: FamilySource,
    draft: PedigreeDraft,
    *,
    user: str,
    comment: str,
    flags: Sequence[str] = (),
    now: Optional[datetime] = None,
) -> EditRecord:
    """Write ``draft`` to every cohort pedigree listing the family, then the rest.

    In order, under one lock per data directory:

    1. refuse: a viewer, an empty comment, no change, an error from
       :func:`problems`, a live ghfc-ngs run of one of the cohorts, or a cohort
       file changed since ``source`` was read;
    2. write every new file beside its target (a permission problem stops the
       save here, before anything changed);
    3. rename them over the cohort pedigrees, then the family pedigree;
    4. keep the rows before and after, and append the log and index lines.
    """
    _check_role()
    return _write_edit(source, draft, user=user, comment=comment, flags=flags, now=now)


async def save_pedigree_edit_async(
    source: FamilySource,
    draft: PedigreeDraft,
    *,
    user: str,
    comment: str,
    flags: Sequence[str] = (),
) -> EditRecord:
    """:func:`save_pedigree_edit` for a page: the role is checked on the event
    loop (it reads the session), the files are written in a thread."""
    _check_role()
    return await asyncio.to_thread(
        _write_edit, source, draft, user=user, comment=comment, flags=flags
    )


def _check_role() -> None:
    # Checked here as well as on the button: an event can be crafted.
    if not can_write():
        raise PermissionError("Only curators and administrators can edit a pedigree")


def _write_edit(
    source: FamilySource,
    draft: PedigreeDraft,
    *,
    user: str,
    comment: str,
    flags: Sequence[str] = (),
    now: Optional[datetime] = None,
) -> EditRecord:
    if not user:
        raise ValueError("No user to record the edit under")
    comment = comment.strip()
    if not comment:
        raise PedigreeEditError("Say why the pedigree is edited: a comment is required")
    original = source.original
    changes = describe_changes(original, draft)
    if not changes:
        raise PedigreeEditError("Nothing changed")
    found = problems(original, draft, source.other_ids)
    if has_errors(found):
        raise PedigreeEditError("; ".join(p.message for p in found if p.level == ERROR))
    rerun = rerun_needed(original, draft)
    now = now or datetime.now()
    family_id = source.family_id
    data_dir = source.data_dir

    lock_path = index_file(data_dir).parent / LOCK_FILENAME
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            for f in source.files:
                if (f.path.parent / PROGRESS_FILENAME).exists():
                    raise PedigreeEditError(
                        f"A ghfc-ngs run of {_cohort_name(f.path)} is in progress: "
                        "save once it has finished"
                    )
                current = f.path.read_text(encoding="utf-8")
                if hashlib.sha256(current.encode("utf-8")).hexdigest() != f.sha256:
                    raise PedigreeEditError(
                        f"{f.path.name} changed since the editor was opened: "
                        "reopen it and make the edit again"
                    )
            new_texts = [render_cohort_text(f, original, draft) for f in source.files]
            family_path = (
                get_family_path(data_dir, family_id) / f"{family_id}.pedigree.tsv"
            )
            history = edits_dir(data_dir, family_id)
            edit_id = _edit_id(history, family_id, now)
            before = family_file_text(source.files[0].text, family_id)
            after = family_file_text(new_texts[0], family_id)
            targets = [
                *(f.path for f in source.files),
                family_path,
                history / f"{family_id}.{edit_id}.before.tsv",
                history / f"{family_id}.{edit_id}.after.tsv",
            ]
            temps: List[Path] = []
            try:
                for target, text in zip(targets, [*new_texts, after, before, after]):
                    temps.append(_write_temp(target, text))
            except OSError as e:
                for t in temps:
                    t.unlink(missing_ok=True)
                raise PedigreeEditError(
                    f"Cannot write {target}: {e.strerror or e}"
                ) from e
            for temp, target in zip(temps, targets):
                os.replace(temp, target)

            record = EditRecord(
                family_id=family_id,
                edit_id=edit_id,
                timestamp=now.isoformat(timespec="seconds"),
                user=user,
                comment=comment,
                changes=tuple(changes),
                cohorts=tuple(source.cohorts),
                rerun=rerun,
                flags=tuple(flags),
            )
            shared = [record.timestamp, user, comment, "; ".join(changes),
                      ",".join(record.cohorts), rerun]  # fmt: skip
            append_tsv_row(
                history / f"{family_id}.pedigree_edits.tsv",
                LOG_HEADER,
                [edit_id, *shared, ",".join(flags)],
            )
            append_tsv_row(
                index_file(data_dir), INDEX_HEADER, [family_id, edit_id, *shared]
            )
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    return record


# ---------------------------------------------------------------------------
# Reading the history
# ---------------------------------------------------------------------------


def _records(path: Path, family_id: Optional[str] = None) -> List[EditRecord]:
    if not path.exists():
        return []
    out: List[EditRecord] = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
            out.append(
                EditRecord(
                    family_id=row.get("FID") or family_id or "",
                    edit_id=row.get("EditId") or "",
                    timestamp=row.get("Timestamp") or "",
                    user=row.get("User") or "",
                    comment=row.get("Comment") or "",
                    changes=tuple(
                        c for c in (row.get("Changes") or "").split("; ") if c
                    ),
                    cohorts=tuple(
                        c for c in (row.get("Cohorts") or "").split(",") if c
                    ),
                    rerun=row.get("Rerun") or "",
                    flags=tuple(c for c in (row.get("Flags") or "").split(",") if c),
                )
            )
    return out


def read_family_edits(data_dir: Path, family_id: str) -> List[EditRecord]:
    """A family's edits, oldest first, from its own log."""
    path = edits_dir(data_dir, family_id) / f"{family_id}.pedigree_edits.tsv"
    return _records(path, family_id)


def last_edits(data_dir: Path) -> Dict[str, EditRecord]:
    """The last edit of each family, from the data directory's index."""
    return {r.family_id: r for r in _records(index_file(data_dir))}


@dataclass(frozen=True)
class EditStatus:
    """A family's last pedigree edit, and whether a ghfc-ngs run followed it."""

    last: EditRecord
    #: An edit that changes what the steps read (not a phenotype alone) was
    #: saved after the last run started, or the cohort has no run record.
    awaiting_rerun: bool


def _as_aware(timestamp: str) -> Optional[datetime]:
    try:
        moment = datetime.fromisoformat(timestamp)
    except ValueError:
        return None
    return moment if moment.tzinfo else moment.astimezone()


def edit_statuses(
    records: Iterable[EditRecord], last_run_started: Optional[datetime]
) -> Dict[str, EditStatus]:
    """Per family: its last edit, and whether one still waits for a run.

    The run record says when the cohort's last ghfc-ngs run started; an edit
    saved after it (or with no record at all) is not in any step's outputs yet.
    """
    out: Dict[str, EditStatus] = {}
    for record in records:
        moment = _as_aware(record.timestamp)
        waits = record.rerun != RERUN_NONE and (
            last_run_started is None or moment is None or moment > last_run_started
        )
        earlier = out.get(record.family_id)
        out[record.family_id] = EditStatus(
            record, waits or (earlier is not None and earlier.awaiting_rerun)
        )
    return out


def read_edit_index(data_dir: Path) -> List[EditRecord]:
    """Every edit of the data directory, oldest first, from its index."""
    return _records(index_file(data_dir))


def draft_from_text(text: str, family_id: str) -> PedigreeDraft:
    """The family's rows in a kept ``before``/``after`` file, as a draft."""
    lines = text.splitlines(keepends=True)
    draft = PedigreeDraft(family_id)
    for index, line in enumerate(lines):
        if line.startswith("#") or not line.strip():
            continue
        cells = _split(line)
        if index == 0 and cells[0].strip() == "FID":
            continue
        if len(cells) >= 6 and cells[0].strip() == family_id:
            cells = [c.strip() for c in cells[:6]] + cells[6:]
            if cells[1] not in draft:
                draft.members.append(
                    Member(
                        cells[1],
                        cells[2],
                        cells[3],
                        cells[4],
                        cells[5],
                        tuple(cells[6:]),
                    )
                )
    return draft


def edit_rows(data_dir: Path, family_id: str, edit_id: str, which: str) -> str:
    """The family's rows ``before`` or ``after`` an edit, as kept."""
    if which not in ("before", "after"):
        raise ValueError(which)
    return (
        edits_dir(data_dir, family_id) / f"{family_id}.{edit_id}.{which}.tsv"
    ).read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Steps still waiting for a rerun
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PendingStep:
    """A step whose last run read other pedigree rows than the current ones."""

    step: str
    label: str
    path: Path


def _snapshot(path: Path) -> Tuple[Dict[str, str], Dict[str, Tuple[str, ...]]]:
    header: Dict[str, str] = {}
    rows: Dict[str, Tuple[str, ...]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            key, _, value = line[1:].partition("\t")
            header[key] = value
        elif line.strip():
            cells = line.split("\t")
            if len(cells) >= 5:
                rows.setdefault(cells[1], tuple(cells[:5]))
    return header, rows


def family_rows(cohort_text: str, family_id: str) -> Dict[str, Tuple[str, ...]]:
    """The family's rows as the steps' snapshots hold them: normalised, first per member."""
    rows: Dict[str, Tuple[str, ...]] = {}
    for index, line in enumerate(cohort_text.splitlines()):
        if line.startswith("#") or not line.strip():
            continue
        cells = line.split("\t")
        if index == 0 and cells[0].strip() == "FID":
            continue
        if len(cells) >= 6 and cells[0].strip() == family_id:
            row = snapshot_row(cells)
            rows.setdefault(row[1], row)
    return rows


#: Steps whose snapshot leaves dummy rows out (main.nf: the STR table never reads them).
_NO_DUMMIES = {"expansionhunter"}


def pending_steps(
    data_dir: Path, family_id: str, cohort_text: str
) -> List[PendingStep]:
    """The family's step outputs built from other pedigree rows than the cohort's now.

    Every ``*.inputs.tsv`` under the family directory is a step's snapshot
    (``#step`` names it). Rows are compared by member, so the Wombat snapshot's
    order (sorted by sample) does not matter. The comparison is the pipeline's
    own: when it differs, the next run of that step rebuilds the family.
    """
    current = family_rows(cohort_text, family_id)
    out: List[PendingStep] = []
    family_dir = get_family_path(data_dir, family_id)
    if not family_dir.is_dir():
        return out
    for path in sorted(family_dir.rglob("*.inputs.tsv")):
        if EDITS_DIRNAME in path.parts:
            continue
        try:
            header, rows = _snapshot(path)
        except (OSError, UnicodeDecodeError):
            continue
        step = header.get("step")
        if not step or header.get("family", family_id) != family_id:
            continue
        expected = current
        if step in _NO_DUMMIES:
            expected = {k: v for k, v in current.items() if not is_dummy_id(k)}
        if rows != expected:
            out.append(PendingStep(step, header.get("label", ""), path))
    return out
