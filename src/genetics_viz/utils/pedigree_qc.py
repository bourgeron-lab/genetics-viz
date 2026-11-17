"""Sample QC against a pedigree draft: what the pipeline would flag, and fixes.

The ghfc-ngs ``sampleqc`` step compares each within-family pair's observed
relationship (KING kinship, IBS0) with the one the pedigree implies, and flags
the contradictions. The editor needs the same verdict for a pedigree that is
not saved yet, so this module mirrors ``modules/sampleqc/scripts/sampleqc_classify``
(ghfc-ngs ``163ce59``): :class:`Kinship` is its ``Pedigree`` (founders
unrelated, a parent with no row a founder, ``0`` nobody), :func:`compare` its
``compare``. The observed side is read as written (``pairs.tsv``), so only the
pedigree side changes.

:func:`suggest` turns the family's QC into candidate edits: the sex of a
karyotype, dummy parents for siblings the pedigree does not join, a dummy for a
parent who is not one, the parent half or full sibs should (not) share, swaps.
Each candidate is applied to a copy of the draft and kept only when it lowers
the contradictions (or, at equal contradictions, the unexplained relatives):
what is offered is what the pipeline would agree with.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .pedigree_edit import (
    FATHER,
    FEMALE,
    MALE,
    MOTHER,
    UNKNOWN,
    PedigreeDraft,
    PedigreeEditError,
    norm_sex,
)
from .sample_qc import FlagRow, PairRow

#: KING's classes, as ``sampleqc_classify`` cuts them.
KING_DUPLICATE = 0.354
KING_FIRST = 0.177
KING_SECOND = 0.0884
KING_THIRD = 0.0442
PO_IBS0_IBS2_MAX = 0.01

DEGREE = {
    "duplicate": 0,
    "parent-offspring": 1, "full-sibs": 1, "first-degree": 1,
    "half-sibs": 2, "second-degree": 2,
    "third-degree": 3,
    "unrelated": 4,
}  # fmt: skip

#: What the pipeline flags as information, not as a contradiction.
INFO_CODES = {"RELATED_BEYOND_PEDIGREE", "CONSANGUINITY", "MZ_TWIN_CANDIDATE"}

#: A karyotype's pedigree sex; an atypical one (X0, XXY…) suggests nothing.
KARYOTYPE_SEX = {"XX": FEMALE, "XY": MALE}

#: How many swaps to offer at most: the best ones.
MAX_SWAPS = 3


def kinship_class(phi: float, ibs0_ibs2: Optional[float]) -> str:
    """``sampleqc_classify.kinship_class``: observed (with IBS0) or expected."""
    if phi > KING_DUPLICATE:
        return "duplicate"
    if phi > KING_FIRST:
        if ibs0_ibs2 is None:
            return "first-degree"
        return "parent-offspring" if ibs0_ibs2 < PO_IBS0_IBS2_MAX else "full-sibs"
    if phi > KING_SECOND:
        return "second-degree"
    if phi > KING_THIRD:
        return "third-degree"
    return "unrelated"


class Kinship:
    """``sampleqc_classify.Pedigree`` over a draft's members and implied dummies."""

    def __init__(self, draft: PedigreeDraft) -> None:
        self.parents: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
        for m in [*draft.members, *draft.implied]:
            self.parents[m.iid] = tuple(
                None if p == UNKNOWN else p
                for p in (m.parent(FATHER), m.parent(MOTHER))
            )  # type: ignore[assignment]
        for father, mother in list(self.parents.values()):
            for p in (father, mother):
                if p is not None and p not in self.parents:
                    self.parents[p] = (None, None)
        self.depth: Dict[str, int] = {}
        for iid in self.parents:
            self._depth(iid, set())
        self.memo: Dict[Tuple[str, str], float] = {}

    def _depth(self, iid: str, seen: Set[str]) -> int:
        if iid in self.depth:
            return self.depth[iid]
        if iid in seen:
            raise PedigreeEditError(f"the pedigree has a cycle through {iid}")
        seen.add(iid)
        parents = [p for p in self.parents[iid] if p is not None]
        d = 1 + max(self._depth(p, seen) for p in parents) if parents else 0
        seen.discard(iid)
        self.depth[iid] = d
        return d

    def kinship(self, a: Optional[str], b: Optional[str]) -> float:
        if a is None or b is None:
            return 0.0
        key = (a, b) if a <= b else (b, a)
        if key in self.memo:
            return self.memo[key]
        if a == b:
            father, mother = self.parents[a]
            value = 0.5 * (1 + self.kinship(father, mother))
        else:
            if self.depth[a] < self.depth[b]:
                a, b = b, a
            father, mother = self.parents[a]
            value = 0.5 * (self.kinship(father, b) + self.kinship(mother, b))
        self.memo[key] = value
        return value

    def relationship(self, a: str, b: str) -> str:
        pa, pb = self.parents.get(a, (None, None)), self.parents.get(b, (None, None))
        if a in pb or b in pa:
            return "parent-offspring"
        shared = {p for p in pa if p is not None} & {p for p in pb if p is not None}
        if len(shared) == 2:
            return "full-sibs"
        if len(shared) == 1:
            return "half-sibs"
        return kinship_class(self.kinship(a, b), None)

    def siblings(self, a: str, b: str) -> bool:
        return bool(
            {p for p in self.parents[a] if p} & {p for p in self.parents[b] if p}
        )

    def may_be_full_sibs(self, a: str, b: str) -> bool:
        (fa, ma), (fb, mb) = self.parents[a], self.parents[b]
        if fa and fa == fb:
            return ma is None or mb is None
        if ma and ma == mb:
            return fa is None or fb is None
        return False

    def both_parents(self, iid: str) -> bool:
        return all(self.parents[iid])

    def spouses(self, a: str, b: str) -> bool:
        return any(set(parents) == {a, b} for parents in self.parents.values())


@dataclass(frozen=True)
class PairCheck:
    """A within-family pair as the pipeline would read it against a draft."""

    a: str
    b: str
    expected: str
    observed: str
    kinship: Optional[float]
    #: The flag the pipeline would raise, ``None`` when the two agree.
    code: Optional[str] = None
    severity: Optional[str] = None

    @property
    def contradicts(self) -> bool:
        return self.code is not None and self.code not in INFO_CODES

    @property
    def unexplained(self) -> bool:
        """Related beyond the pedigree: a parent it lacks may be shared."""
        return self.code == "RELATED_BEYOND_PEDIGREE"

    def describe(self) -> str:
        phi = (
            f"kinship {self.kinship:.3f}" if self.kinship is not None else "no kinship"
        )
        return f"{self.a} – {self.b}: expected {self.expected}, observed {self.observed} ({phi})"


def compare(
    pair: PairRow,
    ped: Kinship,
    sexes: Dict[str, str],
    confirmed: Set[Tuple[str, str]],
) -> Optional[Tuple[str, str]]:
    """``sampleqc_classify.compare``'s flag and severity, ``None`` when they agree.

    ``confirmed`` holds the (parent, child) pairs observed parent-offspring with
    no flag: full sibs confirmed the children of both parents are not
    HALF_SIB_NOT_FULL.
    """
    a, b = pair.sample_a, pair.sample_b
    observed = pair.observed or "low_confidence"
    if observed == "low_confidence":
        return None
    expected = ped.relationship(a, b)
    if observed == "duplicate":
        if ped.siblings(a, b) and sexes.get(a) == sexes.get(b) != UNKNOWN:
            return "MZ_TWIN_CANDIDATE", "info"
        return "DUPLICATE_OR_MZ", "high"
    exp_d, obs_d = DEGREE[expected], DEGREE.get(observed, 4)
    if (
        expected == observed
        or (exp_d == obs_d and exp_d != 1)
        or (expected == "first-degree" and obs_d == 1)
    ):
        return None
    if (
        expected == "half-sibs"
        and observed == "full-sibs"
        and ped.may_be_full_sibs(a, b)
    ):
        return None
    if obs_d < exp_d and not (
        obs_d == 1 and ped.both_parents(a) and ped.both_parents(b)
    ):
        if ped.spouses(a, b):
            return "CONSANGUINITY", "info"
        return "RELATED_BEYOND_PEDIGREE", "info"
    if expected == "parent-offspring":
        return "NOT_PARENT", "high"
    if expected == "full-sibs" and observed == "second-degree":
        parents = [p for p in ped.parents[a] if p]
        if all((p, c) in confirmed for p in parents for c in (a, b)):
            return None
        return "HALF_SIB_NOT_FULL", "high"
    if expected == "half-sibs" and observed == "full-sibs":
        return "FULL_SIB_NOT_HALF", "high"
    if expected == "unrelated":
        return "UNEXPECTED_RELATED_IN_FAMILY", "high" if obs_d <= 2 else "medium"
    if observed == "unrelated":
        return "EXPECTED_RELATED_UNRELATED", "high" if exp_d <= 2 else "medium"
    severity = "high" if min(exp_d, obs_d) == 1 or abs(exp_d - obs_d) > 1 else "medium"
    return "RELATIONSHIP_MISMATCH", severity


def within_family(pairs: Iterable[PairRow], draft: PedigreeDraft) -> List[PairRow]:
    """The pairs of two sequenced members of the draft."""
    members = set(draft.sequenced_ids)
    return [
        p
        for p in pairs
        if not p.cross_family and p.sample_a in members and p.sample_b in members
    ]


def check_pairs(
    draft: PedigreeDraft,
    pairs: Sequence[PairRow],
    upd_pairs: Iterable[Tuple[str, str]] = (),
) -> Optional[List[PairCheck]]:
    """Every within-family pair against the draft; ``None`` for a cyclic draft.

    ``upd_pairs`` (parent, child) are pairs the pipeline flags ``UPD_CANDIDATE``:
    their opposite homozygotes lie on one chromosome, so they are not a
    parent-offspring contradiction.
    """
    try:
        ped = Kinship(draft)
    except PedigreeEditError:
        return None
    sexes = {m.iid: norm_sex(m.sex) for m in draft.members}
    upd = {frozenset(p) for p in upd_pairs}
    confirmed: Set[Tuple[str, str]] = set()
    rows = within_family(pairs, draft)
    for p in rows:
        if p.observed == "parent-offspring":
            for parent, child in ((p.sample_a, p.sample_b), (p.sample_b, p.sample_a)):
                if parent in ped.parents.get(child, ()):
                    if compare(p, ped, sexes, set()) is None:
                        confirmed.add((parent, child))
    out: List[PairCheck] = []
    for p in rows:
        found = compare(p, ped, sexes, confirmed)
        if found is not None and found[0] == "NOT_PARENT":
            if frozenset((p.sample_a, p.sample_b)) in upd:
                found = None
        code, severity = found if found else (None, None)
        out.append(
            PairCheck(
                p.sample_a,
                p.sample_b,
                ped.relationship(p.sample_a, p.sample_b),
                p.observed or "low_confidence",
                p.kinship,
                code,
                severity,
            )
        )
    return out


@dataclass(frozen=True)
class Fit:
    """How well a draft agrees with Sample QC: fewer is better, in this order."""

    contradictions: int
    unexplained: int

    def __lt__(self, other: "Fit") -> bool:
        return (self.contradictions, self.unexplained) < (
            other.contradictions,
            other.unexplained,
        )


def fit(
    draft: PedigreeDraft,
    pairs: Sequence[PairRow],
    upd_pairs: Iterable[Tuple[str, str]] = (),
) -> Optional[Fit]:
    checks = check_pairs(draft, pairs, upd_pairs)
    if checks is None:
        return None
    return Fit(sum(c.contradicts for c in checks), sum(c.unexplained for c in checks))


# ---------------------------------------------------------------------------
# Suggestions
# ---------------------------------------------------------------------------


@dataclass
class Suggestion:
    """A candidate edit, with why it is offered and what it does to the fit."""

    key: str
    title: str
    #: The QC evidence: the pairs or the karyotype it answers.
    why: List[str]
    change: Callable[[PedigreeDraft], object]
    #: The flags it answers, recorded with the edit.
    flags: Tuple[str, ...] = ()
    before: Optional[Fit] = None
    after: Optional[Fit] = None
    #: Something the curator should weigh before applying it.
    caution: str = ""
    kind: str = "pairs"  # "sex" suggestions do not change the pairs' fit

    def apply(self, draft: PedigreeDraft) -> None:
        self.change(draft)


def _try(
    draft: PedigreeDraft, change: Callable[[PedigreeDraft], object]
) -> Optional[PedigreeDraft]:
    trial = draft.copy()
    try:
        change(trial)
    except (PedigreeEditError, ValueError):
        return None
    return trial


def _new_dummy(child: str, role: str) -> Callable[[PedigreeDraft], object]:
    def change(d: PedigreeDraft) -> None:
        d.set_parent(child, role, None)
        d.add_dummy_parents([child], roles=[role])

    return change


_FULL = ("full-sibs", "first-degree")


def _components(
    nodes: Sequence[str], edges: Iterable[Tuple[str, str]]
) -> List[List[str]]:
    """Connected components of ``nodes`` (in their order), size 2 or more."""
    group = {n: n for n in nodes}

    def root(n: str) -> str:
        while group[n] != n:
            group[n] = group[group[n]]
            n = group[n]
        return n

    for a, b in edges:
        if a in group and b in group:
            group[root(a)] = root(b)
    out: Dict[str, List[str]] = {}
    for n in nodes:
        out.setdefault(root(n), []).append(n)
    return [c for c in out.values() if len(c) > 1]


def _sibship_suggestions(
    draft: PedigreeDraft,
    checks: Sequence[PairCheck],
    offer: Callable[..., None],
) -> None:
    """Parents siblings should share, from their observed relationship.

    Full sibs lacking the same parent get one shared dummy (several at once);
    a full sib lacking a parent the other has gets that parent; half sibs who
    share one parent and lack the other get one dummy each. With a parent
    written 0 the pipeline flags none of this: the sibship is only undecided.
    """
    po = {
        frozenset((m.iid, m.parent(r))) for m in draft.members for r in (FATHER, MOTHER)
    }
    full = [
        (c.a, c.b)
        for c in checks
        if c.observed in _FULL and frozenset((c.a, c.b)) not in po
    ]
    for roles in ([FATHER, MOTHER], [FATHER], [MOTHER]):
        # Lacking exactly these: one lacking both is in the first group only.
        lacking = [
            m.iid
            for m in draft.members
            if not m.is_dummy
            and all(m.parent(r) == UNKNOWN for r in roles)
            and all(m.parent(r) != UNKNOWN for r in (FATHER, MOTHER) if r not in roles)
        ]
        for group in _components(lacking, full):
            names = " and ".join(roles) if len(roles) == 1 else "parents"
            offer(
                Suggestion(
                    key=f"shared:{'+'.join(group)}:{'+'.join(roles)}",
                    title=f"Give {', '.join(group)} shared dummy {names} (full sibs)",
                    why=[
                        c.describe()
                        for c in checks
                        if {c.a, c.b} <= set(group) and (c.a, c.b) in full
                    ],
                    change=lambda d, g=group, r=roles: d.add_dummy_parents(g, roles=r),
                    flags=("RELATED_BEYOND_PEDIGREE",),
                ),
                settles=True,
            )
    by_pair = {(c.a, c.b): c for c in checks}
    for a, b in full:
        ma, mb = draft.get(a), draft.get(b)
        if ma is None or mb is None:
            continue
        for r in (FATHER, MOTHER):
            for child, other in ((ma, mb), (mb, ma)):
                parent = other.parent(r)
                if child.parent(r) == UNKNOWN and parent != UNKNOWN:
                    offer(
                        Suggestion(
                            key=f"fullsib:{child.iid}:{r}:{parent}",
                            title=f"Make {parent} {child.iid}'s {r} too (full sibs)",
                            why=[by_pair[(a, b)].describe()],
                            change=lambda d, x=child.iid, rr=r, p=parent: d.set_parent(
                                x, rr, p
                            ),
                        ),
                        settles=True,
                    )
    for c in checks:
        if c.observed != "second-degree" or c.expected != "half-sibs":
            continue
        ma, mb = draft.get(c.a), draft.get(c.b)
        if ma is None or mb is None:
            continue
        for r in (FATHER, MOTHER):
            # One has this parent, the other lacks it: a dummy of their own.
            for lacker, other in ((ma, mb), (mb, ma)):
                if lacker.parent(r) == UNKNOWN and other.parent(r) != UNKNOWN:
                    offer(
                        Suggestion(
                            key=f"own:{lacker.iid}:{r}",
                            title=f"Give {lacker.iid} a dummy {r} of their own"
                            f" (half sib of {other.iid})",
                            why=[c.describe()],
                            change=_new_dummy(lacker.iid, r),
                        ),
                        settles=True,
                    )
            if ma.parent(r) == UNKNOWN and mb.parent(r) == UNKNOWN:

                def each(
                    d: PedigreeDraft, x: str = c.a, y: str = c.b, rr: str = r
                ) -> None:
                    d.add_dummy_parents([x], roles=[rr])
                    d.add_dummy_parents([y], roles=[rr])

                offer(
                    Suggestion(
                        key=f"halfeach:{c.a}:{c.b}:{r}",
                        title=f"Give {c.a} and {c.b} a dummy {r} each (half sibs)",
                        why=[c.describe()],
                        change=each,
                    ),
                    settles=True,
                )


def suggest(
    draft: PedigreeDraft,
    pairs: Sequence[PairRow],
    samples: Sequence[Dict[str, object]] = (),
    flags: Sequence[FlagRow] = (),
) -> List[Suggestion]:
    """Edits Sample QC supports for this draft, best first.

    ``pairs`` and ``samples`` are the family's ``pairs.tsv`` and
    ``samples.tsv`` rows; ``flags`` its flags (for ``UPD_CANDIDATE`` pairs).
    A pair suggestion is kept only if it improves :func:`fit`; a sex one when
    the karyotype is XX or XY and the draft says otherwise.
    """
    upd = [
        (f.sample_a, f.sample_b or "")
        for f in flags
        if f.code == "UPD_CANDIDATE" and f.sample_b
    ]
    checks = check_pairs(draft, pairs, upd)
    if checks is None:
        return []
    base = fit(draft, pairs, upd)
    out: List[Suggestion] = []
    seen: Set[str] = set()

    def offer(s: Suggestion, settles: bool = False) -> None:
        """Keep ``s`` if it improves the fit; one that ``settles`` an undecided
        sibship (no flag either way) if it does not make the fit worse."""
        if s.key in seen:
            return
        trial = _try(draft, s.change)
        if trial is None:
            return
        after = fit(trial, pairs, upd)
        if after is None or (base < after if settles else not after < base):  # type: ignore[operator]
            return
        s.before, s.after = base, after
        seen.add(s.key)
        out.append(s)

    contradicted = {m for c in checks if c.contradicts for m in (c.a, c.b)}
    _sibship_suggestions(draft, checks, offer)

    # Sex from the karyotype.
    sex_fixes: List[Suggestion] = []
    for row in samples:
        iid = str(row.get("sample") or "")
        member = draft.get(iid)
        wanted = KARYOTYPE_SEX.get(str(row.get("karyotype") or ""))
        if member is None or wanted is None or norm_sex(member.sex) == wanted:
            continue
        label = "male" if wanted == MALE else "female"
        caution = (
            f"{iid}'s relationships contradict the pedigree too: a sample swap"
            " is more likely than a sex entry error."
            if iid in contradicted
            else ""
        )
        sex_fixes.append(
            Suggestion(
                key=f"sex:{iid}",
                title=f"Record {iid} as {label}",
                why=[
                    f"{iid}: karyotype {row.get('karyotype')} ({row.get('sex_evidence') or 'copy numbers'})"
                ],
                change=lambda d, i=iid, s=wanted: d.set_sex(i, s),
                flags=("SEX_MISMATCH",)
                if norm_sex(member.sex) != UNKNOWN
                else ("SEX_UNKNOWN_IN_PED",),
                caution=caution,
                kind="sex",
            )
        )

    for c in checks:
        if c.code is None:
            continue
        a, b = c.a, c.b
        ma, mb = draft.get(a), draft.get(b)
        if ma is None or mb is None:
            continue
        why = [c.describe()]
        flag = (c.code,)
        if c.code == "RELATED_BEYOND_PEDIGREE" and c.observed == "second-degree":
            for r in (FATHER, MOTHER):
                if ma.parent(r) == UNKNOWN and mb.parent(r) == UNKNOWN:
                    offer(Suggestion(
                        key=f"halfshared:{a}:{b}:{r}",
                        title=f"Give {a} and {b} a shared dummy {r} (half sibs)",
                        why=why, flags=flag,
                        change=lambda d, x=a, y=b, rr=r: d.add_dummy_parents([x, y], roles=[rr]),
                    ))  # fmt: skip
        if c.code == "NOT_PARENT":
            for parent, child in ((a, b), (b, a)):
                member = draft.get(child)
                for r in (FATHER, MOTHER):
                    if member is not None and member.parent(r) == parent:
                        offer(Suggestion(
                            key=f"notparent:{child}:{r}",
                            title=f"Replace {child}'s {r} {parent} with a dummy",
                            why=why, flags=flag,
                            change=_new_dummy(child, r),
                            caution="Non-paternity or a sample swap: check a swap first.",
                        ))  # fmt: skip
        if c.code == "HALF_SIB_NOT_FULL":
            for child in (a, b):
                for r in (FATHER, MOTHER):
                    offer(Suggestion(
                        key=f"halfsib:{child}:{r}",
                        title=f"Give {child} a dummy {r} of their own (half sibs)",
                        why=why, flags=flag,
                        change=_new_dummy(child, r),
                    ))  # fmt: skip
        if c.code == "FULL_SIB_NOT_HALF":
            for r in (FATHER, MOTHER):
                pa, pb = ma.parent(r), mb.parent(r)
                if pa == pb:
                    continue
                for child, parent in ((b, pa), (a, pb)):
                    if parent == UNKNOWN:
                        continue
                    offer(Suggestion(
                        key=f"fullnothalf:{child}:{r}:{parent}",
                        title=f"Make {parent} {child}'s {r} too (full sibs)",
                        why=why, flags=flag,
                        change=lambda d, x=child, rr=r, p=parent: d.set_parent(x, rr, p),
                    ))  # fmt: skip

    # Swaps, when something contradicts the pedigree: the best few.
    if base is not None and base.contradictions:
        swaps: List[Suggestion] = []
        ids = draft.sequenced_ids
        for i, x in enumerate(ids):
            for y in ids[i + 1 :]:
                trial = _try(draft, lambda d, p=x, q=y: d.swap(p, q))
                if trial is None:
                    continue
                after = fit(trial, pairs, upd)
                if after is None or not after < base:
                    continue
                swaps.append(Suggestion(
                    key=f"swap:{x}:{y}", title=f"Swap {x} and {y} (sample swap)",
                    why=[c.describe() for c in checks if c.contradicts and {c.a, c.b} & {x, y}],
                    change=lambda d, p=x, q=y: d.swap(p, q),
                    flags=tuple(sorted({c.code for c in checks if c.contradicts and {c.a, c.b} & {x, y}})),  # type: ignore[misc]
                    before=base, after=after,
                ))  # fmt: skip
        swaps.sort(key=lambda s: (s.after.contradictions, s.after.unexplained))  # type: ignore[union-attr]
        for s in swaps[:MAX_SWAPS]:
            if s.key not in seen:
                seen.add(s.key)
                out.append(s)

    out.sort(key=lambda s: (s.after.contradictions, s.after.unexplained))  # type: ignore[union-attr]
    return sex_fixes + out


@dataclass
class QcContext:
    """What the editor gets from the Sample QC tab, for the fit and suggestions."""

    pairs: List[PairRow] = field(default_factory=list)
    samples: List[Dict[str, object]] = field(default_factory=list)
    flags: List[FlagRow] = field(default_factory=list)

    @property
    def upd_pairs(self) -> List[Tuple[str, str]]:
        return [
            (f.sample_a, f.sample_b or "")
            for f in self.flags
            if f.code == "UPD_CANDIDATE" and f.sample_b
        ]

    def fit(self, draft: PedigreeDraft) -> Optional[Fit]:
        return fit(draft, self.pairs, self.upd_pairs) if self.pairs else None

    def checks(self, draft: PedigreeDraft) -> Optional[List[PairCheck]]:
        return check_pairs(draft, self.pairs, self.upd_pairs) if self.pairs else None

    def suggest(self, draft: PedigreeDraft) -> List[Suggestion]:
        return suggest(draft, self.pairs, self.samples, self.flags)
