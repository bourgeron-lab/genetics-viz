"""Merging families into one, under a new family ID or one of theirs.

Some families are one: a pedigree split by mistake, or one where every member
was entered under a family ID of its own. A merge rewrites every member's
line of the cohort pedigrees under the target ID, and writes the target's
family pedigree as ``FAMILIAL_PEDIGREE`` would. For ghfc-ngs:

- a **new** family ID is a new family: a normal run builds it from the
  existing gVCFs (the stub test ``addfamily``);
- an **existing** one gains members: ``--clean-stale-families``;
- the other families leave the pedigree. Their directories stay on disk,
  orphaned: the cohort tables are rebuilt without them, and the Workflow tab
  lists them for removal by hand.

Everything the app's curators wrote is keyed by family ID: validations,
diagnostics, notes, Sample QC acknowledgements, report conclusions. A merge
**re-keys** their rows to the target ID (:func:`rekey_curation`), after
copying each file it changes into the target's history. A cross-family
acknowledgement (``"FA,FB"``) is left as it is: once the families are one,
the pipeline flags that pair under another code.

The plan is the vault note ``Pedigree editor`` (phase 5).
"""

from __future__ import annotations

import asyncio
import fcntl
import hashlib
import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set

from .diagnostic_files import diagnostic_paths
from .pedigree_edit import (
    ERROR,
    INDEX_HEADER,
    LOCK_FILENAME,
    LOG_HEADER,
    RERUN_CLEAN,
    RERUN_NORMAL,
    EditRecord,
    FamilySource,
    Member,
    PedigreeDraft,
    PedigreeEditError,
    Problem,
    _check_role,
    _cohort_name,
    _edit_id,
    _ending,
    _member_line,
    _split,
    _write_temp,
    describe_changes,
    edits_dir,
    family_file_text,
    index_file,
    live_run,
    outside_ids,
    problems,
)
from .sharding import get_family_path
from .tsv import append_tsv_row, rewrite_tsv_rows
from .validation_files import validation_paths

#: The joint of a suggested family ID: the merged IDs in order.
SEPARATOR = "_"
#: An old family's ``Rerun``: it left the pedigree, merged into another.
RERUN_MERGED = "merged"
#: A family ID: no whitespace, no path separator, no comma (the acks' pairs).
_FID = re.compile(r"^[^\s/\\,]+$")


def suggested_fid(family_ids: Sequence[str]) -> str:
    """The merged families' IDs joined, in the order given."""
    return SEPARATOR.join(family_ids)


def curation_files(data_dir: Path) -> List[Path]:
    """Every curation file whose rows carry a family ID (``FID``)."""
    data_dir = Path(data_dir)
    return [
        *(p for _, p in validation_paths(data_dir)),
        *(p for _, p in diagnostic_paths(data_dir)),
        data_dir / "validations" / "sampleqc.tsv",
        data_dir / "notes" / "notes.tsv",
        data_dir / "reports" / "family_summaries.tsv",
        data_dir / "reports" / "family_status.tsv",
    ]


def curation_counts(data_dir: Path, family_ids: Iterable[str]) -> Dict[Path, int]:
    """Rows of each curation file keyed by one of ``family_ids`` (what a merge moves)."""
    wanted = set(family_ids)
    out: Dict[Path, int] = {}
    for path in curation_files(data_dir):
        if not path.exists():
            continue
        with open(path, encoding="utf-8") as f:
            header = f.readline().rstrip("\r\n").split("\t")
            if "FID" not in header:
                continue
            col = header.index("FID")
            n = sum(
                1
                for line in f
                if (cells := line.rstrip("\r\n").split("\t"))
                and len(cells) > col
                and cells[col] in wanted
            )
        if n:
            out[path] = n
    return out


def rekey_curation(data_dir: Path, old: Iterable[str], target: str) -> Dict[Path, int]:
    """Give every curation row of an ``old`` family ID the ``target`` one.

    Each file is rewritten under its own lock (:func:`tsv.rewrite_tsv_rows`).
    Returns the rows changed per file.
    """
    old_ids = set(old) - {target}
    out: Dict[Path, int] = {}
    for path in curation_files(data_dir):
        count = 0

        def change(rows: List[Dict[str, str]]) -> bool:
            nonlocal count
            for row in rows:
                if row.get("FID") in old_ids:
                    row["FID"] = target
                    count += 1
            return count > 0

        if path.exists() and rewrite_tsv_rows(path, change):
            out[path] = count
    return out


@dataclass
class MergeSource:
    """The families to merge, as read, and the ID they will share."""

    sources: List[FamilySource]
    target: str

    @property
    def data_dir(self) -> Path:
        return self.sources[0].data_dir

    @property
    def family_ids(self) -> List[str]:
        return [s.family_id for s in self.sources]

    @property
    def cohorts(self) -> List[str]:
        return self.sources[0].cohorts

    def original(self) -> PedigreeDraft:
        """Every family's members as one draft, under the target ID."""
        return PedigreeDraft(
            self.target, [m for s in self.sources for m in s.original.members]
        )

    def as_source(self) -> FamilySource:
        """The merge as one family, for the editor: its files and original."""
        return FamilySource(
            self.data_dir, self.target, self.sources[0].files, self.original()
        )

    def other_ids(self) -> Set[str]:
        return outside_ids(self.sources)


def load_merge(
    data_dir: Path, family_ids: Sequence[str], cohort_files: Sequence[Sequence[Path]]
) -> MergeSource:
    """Read the families to merge (each from its cohort pedigrees, which must
    be the same files for all), with the suggested target ID."""
    from .pedigree_edit import load_family

    if len(family_ids) < 2:
        raise PedigreeEditError("Choose two families or more to merge")
    if any(list(f) != list(cohort_files[0]) for f in cohort_files[1:]):
        raise PedigreeEditError(
            "These families are not in the same cohort pedigrees: they cannot"
            " be merged in one edit"
        )
    sources = [
        load_family(data_dir, fid, files)
        for fid, files in zip(family_ids, cohort_files)
    ]
    return MergeSource(sources, suggested_fid(family_ids))


def target_problems(merge: MergeSource, target: str) -> List[Problem]:
    """What is wrong with ``target`` as the merged family's ID."""
    out: List[Problem] = []
    if not target or not _FID.match(target):
        out.append(Problem(ERROR, "The family ID must have no space, comma or slash"))
        return out
    if target in merge.family_ids:
        return out
    taken = any(
        _split(line)[0].strip() == target
        for f in merge.sources[0].files
        for line in f.lines
        if line.strip() and not line.startswith("#")
    )
    if taken:
        out.append(
            Problem(ERROR, f"Family {target} already exists in the cohort pedigree")
        )
    elif get_family_path(merge.data_dir, target).exists():
        out.append(
            Problem(
                ERROR,
                f"A directory for family {target} already exists: choose another ID"
                " (or remove its leftovers first)",
            )
        )
    return out


def merge_problems(
    merge: MergeSource, draft: PedigreeDraft, target: str
) -> List[Problem]:
    """The target's problems, then the merged draft's as one family."""
    original = merge.original()
    return target_problems(merge, target) + problems(original, draft, merge.other_ids())


def merge_changes(merge: MergeSource, draft: PedigreeDraft, target: str) -> List[str]:
    """The merge first, then each member's change (as the pipeline words them)."""
    others = [f for f in merge.family_ids if f != target]
    head = (
        f"{', '.join(others)} merged into {target}"
        if target in merge.family_ids
        else f"{', '.join(merge.family_ids)} merged into {target} (a new family)"
    )
    return [head, *describe_changes(merge.original(), draft)]


def merge_rerun(merge: MergeSource, target: str) -> str:
    """A new family is built by a normal run; an existing one gains members."""
    return RERUN_CLEAN if target in merge.family_ids else RERUN_NORMAL


def render_merge_text(
    text: str, path: Path, merge: MergeSource, draft: PedigreeDraft, target: str
) -> str:
    """The cohort file with the merged families' lines replaced, where the first
    of them stood, by the merged family's; every other line as it was.

    A member unchanged and already under the target ID keeps its line.
    """
    fids = set(merge.family_ids)
    lines = text.splitlines(keepends=True)
    endings = [_ending(line) for line in lines if _ending(line)]
    newline = endings[0] if endings else "\n"
    kept: Dict[str, str] = {}
    for index, line in enumerate(lines):
        if line.startswith("#") or not line.strip():
            continue
        cells = [c.strip() for c in _split(line)]
        if index == 0 and cells[0] == "FID":
            continue
        if cells[0] == target and len(cells) > 1:
            kept.setdefault(cells[1], line if _ending(line) else line + newline)
    before: Dict[str, Member] = {}
    for s in merge.sources:
        if s.family_id == target:
            before.update({m.iid: m for m in s.original.members})
    merged = [
        kept[m.iid] if m == before.get(m.iid) and m.iid in kept
        else _member_line(target, m, newline)
        for m in draft.members
    ]  # fmt: skip
    out: List[str] = []
    placed = False
    for index, line in enumerate(lines):
        cells = [c.strip() for c in _split(line)]
        is_row = line.strip() and not line.startswith("#")
        if is_row and not (index == 0 and cells[0] == "FID") and cells[0] in fids:
            if not placed:
                if out and not _ending(out[-1]):
                    out[-1] += newline
                out.extend(merged)
                placed = True
            continue
        out.append(line)
    return "".join(out)


def save_merge(
    merge: MergeSource,
    draft: PedigreeDraft,
    target: str,
    *,
    user: str,
    comment: str,
    flags: Sequence[str] = (),
    now: Optional[datetime] = None,
) -> List[EditRecord]:
    """Merge the families under ``target`` and re-key their curation.

    Checked and written as :func:`pedigree_edit.save_pedigree_edit` is (role,
    comment, problems, live run, unchanged files, temp files then renames,
    under the data directory's lock). The target's history keeps every old
    family's rows before and the merged rows after, and a copy of each
    curation file the merge rewrites; each old family's history says where it
    went. One record per family, the target's last.
    """
    _check_role()
    return _write_merge(
        merge, draft, target, user=user, comment=comment, flags=flags, now=now
    )


async def save_merge_async(
    merge: MergeSource,
    draft: PedigreeDraft,
    target: str,
    *,
    user: str,
    comment: str,
    flags: Sequence[str] = (),
) -> List[EditRecord]:
    """:func:`save_merge` for a page (role on the loop, files in a thread)."""
    _check_role()
    return await asyncio.to_thread(
        _write_merge, merge, draft, target, user=user, comment=comment, flags=flags
    )


def _write_merge(
    merge: MergeSource,
    draft: PedigreeDraft,
    target: str,
    *,
    user: str,
    comment: str,
    flags: Sequence[str] = (),
    now: Optional[datetime] = None,
) -> List[EditRecord]:
    if not user:
        raise ValueError("No user to record the edit under")
    comment = comment.strip()
    if not comment:
        raise PedigreeEditError(
            "Say why the families are merged: a comment is required"
        )
    target = target.strip()
    draft = PedigreeDraft(target, list(draft.members))
    found = merge_problems(merge, draft, target)
    if any(p.level == ERROR for p in found):
        raise PedigreeEditError("; ".join(p.message for p in found if p.level == ERROR))
    now = now or datetime.now()
    data_dir = merge.data_dir
    paths = [f.path for f in merge.sources[0].files]
    rerun = merge_rerun(merge, target)
    changes = merge_changes(merge, draft, target)

    lock_path = index_file(data_dir).parent / LOCK_FILENAME
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            for s in merge.sources:
                for f in s.files:
                    if running := live_run(f.path.parent):
                        raise PedigreeEditError(
                            f"A ghfc-ngs run of {_cohort_name(f.path)} is in"
                            f" progress ({running}): save once it has finished"
                        )
                    current = f.path.read_text(encoding="utf-8")
                    if hashlib.sha256(current.encode("utf-8")).hexdigest() != f.sha256:
                        raise PedigreeEditError(
                            f"{f.path.name} changed since the editor was opened: "
                            "reopen it and make the edit again"
                        )
            texts_now = [merge.sources[0].files[i].text for i in range(len(paths))]
            new_texts = [
                render_merge_text(text, path, merge, draft, target)
                for text, path in zip(texts_now, paths)
            ]
            fids = merge.family_ids
            edit_id = max(
                _edit_id(edits_dir(data_dir, fid), fid, now) for fid in {*fids, target}
            )
            merged_rows = family_file_text(new_texts[0], target)
            befores = {fid: family_file_text(texts_now[0], fid) for fid in fids}
            header = (
                merged_rows.splitlines(keepends=True)[:1]
                if merged_rows.startswith("FID")
                else []
            )
            all_before = "".join(header) + "".join(
                "".join(t.splitlines(keepends=True)[len(header) :])
                for t in befores.values()
            )
            history = edits_dir(data_dir, target)
            targets: List[Path] = [*paths]
            texts: List[str] = [*new_texts]
            targets += [
                get_family_path(data_dir, target) / f"{target}.pedigree.tsv",
                history / f"{target}.{edit_id}.before.tsv",
                history / f"{target}.{edit_id}.after.tsv",
            ]
            texts += [merged_rows, all_before, merged_rows]
            for fid in fids:
                if fid == target:
                    continue
                old = edits_dir(data_dir, fid)
                targets += [
                    old / f"{fid}.{edit_id}.before.tsv",
                    old / f"{fid}.{edit_id}.after.tsv",
                ]
                texts += [befores[fid], merged_rows]
            # The curation files the merge rewrites, copied first.
            backup = history / f"{target}.{edit_id}.curation"
            to_rekey = [f for f in fids if f != target]
            for path in curation_counts(data_dir, to_rekey):
                targets.append(backup / path.relative_to(data_dir))
                texts.append(path.read_text(encoding="utf-8"))
            temps: List[Path] = []
            try:
                for tgt, text in zip(targets, texts):
                    temps.append(_write_temp(tgt, text))
            except OSError as e:
                for t in temps:
                    t.unlink(missing_ok=True)
                raise PedigreeEditError(f"Cannot write {tgt}: {e.strerror or e}") from e
            for temp, tgt in zip(temps, targets):
                os.replace(temp, tgt)
            moved = rekey_curation(data_dir, to_rekey, target)

            timestamp = now.isoformat(timespec="seconds")
            records: List[EditRecord] = []
            moved_text = ", ".join(
                f"{n} in {p.relative_to(data_dir)}" for p, n in moved.items()
            )
            for fid in [f for f in fids if f != target] + [target]:
                if fid == target:
                    clauses = [*changes]
                    if moved_text:
                        clauses.append(f"curation re-keyed: {moved_text}")
                    fid_rerun = rerun
                else:
                    clauses = [f"merged into {target}"]
                    fid_rerun = RERUN_MERGED
                record = EditRecord(
                    family_id=fid,
                    edit_id=edit_id,
                    timestamp=timestamp,
                    user=user,
                    comment=comment,
                    changes=tuple(clauses),
                    cohorts=tuple(merge.cohorts),
                    rerun=fid_rerun,
                    flags=tuple(flags),
                )
                shared = [timestamp, user, comment, "; ".join(clauses),
                          ",".join(record.cohorts), fid_rerun]  # fmt: skip
                append_tsv_row(
                    edits_dir(data_dir, fid) / f"{fid}.pedigree_edits.tsv",
                    LOG_HEADER,
                    [edit_id, *shared, ",".join(flags)],
                )
                append_tsv_row(
                    index_file(data_dir), INDEX_HEADER, [fid, edit_id, *shared]
                )
                records.append(record)
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
    return records


@dataclass(frozen=True)
class Orphan:
    """A family merged into another whose directory is still on disk."""

    family_id: str
    into: str
    timestamp: str
    path: Path


def orphaned_families(
    data_dir: Path, records: Iterable[EditRecord], in_pedigree: Iterable[str]
) -> List[Orphan]:
    """Families a merge took out of the pedigree whose directory remains.

    ghfc-ngs leaves them: once the merged family is built, they can go.
    """
    present = set(in_pedigree)
    out: Dict[str, Orphan] = {}
    for r in records:
        if r.rerun != RERUN_MERGED or r.family_id in present:
            continue
        path = get_family_path(data_dir, r.family_id)
        if path.exists():
            into = r.changes[0].removeprefix("merged into ") if r.changes else ""
            out[r.family_id] = Orphan(r.family_id, into, r.timestamp, path)
    return sorted(out.values(), key=lambda o: o.family_id)
