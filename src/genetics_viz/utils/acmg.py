"""ACMG/AMP criteria, points and class (``config/acmg.yaml``).

A diagnostic may record the criteria the curator applied (``Criteria``, e.g.
``PVS1;PS2;PM2_Supporting``), their points (``Points``) and the rules used
(``Ruleset``). The class is computed from the points (Tavtigian et al. 2020),
the curator can choose another: :func:`override` tells it from the record, so
no column says so.

A criterion is written as its code at its own strength, else with the strength
as a suffix (``PM2_Supporting``, ``PVS1_Strong``). An origin may follow,
``:computed`` (by the pipeline, once pywombat writes criteria) or ``:curator``;
none means the curator. Records without criteria are ``legacy``: the verdict
alone, as before v0.46.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

import yaml

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "acmg.yaml"
with open(_CONFIG_PATH) as _f:
    _CONFIG: Dict[str, Any] = yaml.safe_load(_f) or {}

RULESET: str = str(_CONFIG.get("ruleset", ""))
#: Strength → points (benign criteria count them negative).
STRENGTHS: Dict[str, int] = {
    str(k): int(v) for k, v in (_CONFIG.get("strengths") or {}).items()
}
STAND_ALONE_STRENGTH = "StandAlone"
STAND_ALONE: frozenset = frozenset(_CONFIG.get("stand_alone") or ())
#: ``(verdict, minimum points)``, strongest first.
CLASSES: List[tuple] = [
    (str(c["verdict"]), int(c["min"])) for c in _CONFIG.get("classes") or []
]
ORIGINS = ("curator", "computed")
#: Strength names as a curator reads them.
STRENGTH_LABELS = {
    "Supporting": "Supporting",
    "Moderate": "Moderate",
    "Strong": "Strong",
    "VeryStrong": "Very strong",
    STAND_ALONE_STRENGTH: "Stand-alone",
}


@dataclass(frozen=True)
class Criterion:
    code: str
    direction: str  # "pathogenic" or "benign"
    strength: str  # its own strength (Richards 2015)
    text: str
    default: str  # the strength proposed when ticked

    @property
    def strengths(self) -> List[str]:
        """The strengths it may be applied at."""
        if self.strength == STAND_ALONE_STRENGTH:
            return [STAND_ALONE_STRENGTH]
        return list(STRENGTHS)


CRITERIA: Dict[str, Criterion] = {
    code: Criterion(
        code=code,
        direction=str(c["direction"]),
        strength=str(c["strength"]),
        text=str(c.get("text", "")),
        default=str(c.get("default") or c["strength"]),
    )
    for code, c in (_CONFIG.get("criteria") or {}).items()
}


@dataclass(frozen=True)
class Applied:
    """One criterion as applied to a variant."""

    code: str
    strength: str
    origin: str = "curator"

    @property
    def criterion(self) -> Criterion:
        return CRITERIA[self.code]

    @property
    def name(self) -> str:
        """``PM2_Supporting``, or ``PM2`` at its own strength."""
        if self.strength == self.criterion.strength:
            return self.code
        return f"{self.code}_{self.strength}"

    @property
    def points(self) -> int:
        if self.strength == STAND_ALONE_STRENGTH:
            return 0
        sign = -1 if self.criterion.direction == "benign" else 1
        return sign * STRENGTHS[self.strength]

    def __str__(self) -> str:
        return self.name if self.origin == "curator" else f"{self.name}:{self.origin}"


def apply(
    code: str, strength: Optional[str] = None, origin: str = "curator"
) -> Applied:
    """A criterion at ``strength`` (its proposed one by default), checked."""
    if code not in CRITERIA:
        raise ValueError(f"Unknown ACMG criterion: {code!r}")
    criterion = CRITERIA[code]
    strength = strength or criterion.default
    if strength not in criterion.strengths:
        raise ValueError(f"{code} cannot be applied at strength {strength!r}")
    if origin not in ORIGINS:
        raise ValueError(f"Unknown criterion origin: {origin!r}")
    return Applied(code, strength, origin)


def parse_one(text: str) -> Applied:
    """``PM2_Supporting:computed`` → :class:`Applied`."""
    name, _, origin = text.strip().partition(":")
    code, _, strength = name.partition("_")
    if code not in CRITERIA:
        raise ValueError(f"Unknown ACMG criterion: {text!r}")
    return apply(code, strength or CRITERIA[code].strength, origin or "curator")


def parse_criteria(text: Optional[str]) -> List[Applied]:
    """A ``Criteria`` cell → the criteria, in :func:`ordered` order."""
    parts = [p for p in (text or "").split(";") if p.strip()]
    return ordered(parse_one(p) for p in parts)


def ordered(applied: Iterable[Applied]) -> List[Applied]:
    """The config's order (pathogenic strongest first), one per code.

    A code given twice keeps its last strength: the dialog's state is a dict
    keyed by code, but a hand-written cell could repeat one.
    """
    by_code = {a.code: a for a in applied}
    return [by_code[c] for c in CRITERIA if c in by_code]


def format_criteria(applied: Iterable[Applied]) -> str:
    """The ``Criteria`` cell: ``PVS1;PS2;PM2_Supporting``."""
    return ";".join(str(a) for a in ordered(applied))


def total_points(applied: Iterable[Applied]) -> int:
    return sum(a.points for a in applied)


def classify(applied: Sequence[Applied]) -> str:
    """The class the criteria give (``""`` for none).

    A stand-alone benign criterion (BA1) makes the variant benign whatever
    else applies, as in Richards 2015 and the ClinGen BA1 rule.
    """
    if not applied:
        return ""
    if any(a.code in STAND_ALONE for a in applied):
        return CLASSES[-1][0] if CLASSES else ""
    points = total_points(applied)
    for verdict, minimum in CLASSES:
        if points >= minimum:
            return verdict
    return CLASSES[-1][0] if CLASSES else ""


@dataclass(frozen=True)
class Evidence:
    """What a diagnostic record says of its criteria."""

    criteria: List[Applied]
    points: Optional[int]
    ruleset: str
    computed_class: str

    @property
    def recorded(self) -> bool:
        return bool(self.criteria)


def record_evidence(record: Mapping[str, Any]) -> Evidence:
    """Read ``Criteria``, ``Points`` and ``Ruleset`` from a diagnostic record.

    An unreadable cell (a criterion this version does not know) gives no
    criteria rather than an error: the verdict still stands.
    """
    try:
        criteria = parse_criteria(record.get("Criteria"))
    except ValueError:
        criteria = []
    try:
        points: Optional[int] = int(str(record.get("Points") or "").strip())
    except ValueError:
        points = total_points(criteria) if criteria else None
    return Evidence(
        criteria=criteria,
        points=points,
        ruleset=str(record.get("Ruleset") or "").strip(),
        computed_class=classify(criteria),
    )


def override(record: Mapping[str, Any]) -> str:
    """The class the criteria gave when the curator recorded another, else ``""``."""
    evidence = record_evidence(record)
    if not evidence.recorded:
        return ""
    verdict = str(record.get("Diagnostic") or "")
    return evidence.computed_class if evidence.computed_class != verdict else ""


def evidence_columns(applied: Sequence[Applied]) -> List[str]:
    """``[Criteria, Points, Ruleset]`` for a new record (empty without criteria)."""
    applied = ordered(applied)
    if not applied:
        return ["", "", ""]
    return [format_criteria(applied), str(total_points(applied)), RULESET]


def summary(record: Mapping[str, Any], label=None) -> str:
    """``PVS1, PS2, PM2_Supporting · 13 points`` (``""`` without criteria).

    ``label`` words a verdict for the override note (``verdict_label``).
    """
    evidence = record_evidence(record)
    if not evidence.recorded:
        return ""
    text = ", ".join(a.name for a in evidence.criteria)
    if evidence.points is not None:
        unit = "point" if abs(evidence.points) == 1 else "points"
        text += f" · {evidence.points} {unit}"
    if computed := override(record):
        text += f" (criteria give {(label or str)(computed).lower()})"
    return text
