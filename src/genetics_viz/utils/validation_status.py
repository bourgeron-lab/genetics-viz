"""Validation statuses: one row's consensus, and the status filter of the tables.

A row (a variant in one sample) can carry several curators' validations. They
are folded into one status by :func:`consensus`; a row nobody validated has an
empty ``Validation`` and is shown and filtered as :data:`UNVALIDATED`.

The Search page, the Wombat and SV tabs and the per-file validation page all
filter on these statuses, so a new status is an edit here.
"""

from typing import Any, Collection, Dict, Iterable, List, Tuple

#: Filter label and chart category of a row with no (non-ignored) validation.
UNVALIDATED = "TODO"
#: Curators disagree on the row.
CONFLICTING = "conflicting"
#: An SNV validated present as part of a multi-nucleotide variant in phase.
IN_PHASE_MNV = "in phase MNV"

#: Statuses offered by the validation filter menus, in menu order.
FILTER_STATUSES: List[str] = [
    "present",
    "absent",
    "uncertain",
    CONFLICTING,
    UNVALIDATED,
]
#: The same, for pages that list validated rows only.
VALIDATED_STATUSES: List[str] = [s for s in FILTER_STATUSES if s != UNVALIDATED]

#: Menu icon and colour (Quasar names) per filter status.
STATUS_ICONS: Dict[str, Tuple[str, str]] = {
    "present": ("check_circle", "green"),
    "absent": ("cancel", "red"),
    "uncertain": ("help", "orange"),
    CONFLICTING: ("bolt", "amber-9"),
    UNVALIDATED: ("assignment", "grey"),
}

#: Chart colour per filter status (clear of the greys of cytoband stains).
STATUS_COLORS: Dict[str, str] = {
    "present": "#22c55e",
    "absent": "#ef4444",
    "uncertain": "#f59e0b",
    CONFLICTING: "#fb923c",
    UNVALIDATED: "#8b5cf6",
}

_PRESENT = ("present", IN_PHASE_MNV)


def row_status(row: Dict[str, Any]) -> str:
    """The row's status as filters and charts count it."""
    return row.get("Validation") or UNVALIDATED


def status_selected(status: str | None, selected: Collection[str]) -> bool:
    """Whether a row of this status passes a filter selection.

    An empty selection filters nothing. An in-phase MNV counts as ``present``:
    the menu has no entry of its own for it.
    """
    if not selected:
        return True
    status = status or UNVALIDATED
    return status in selected or (status == IN_PHASE_MNV and "present" in selected)


def filter_by_status(
    rows: Iterable[Dict[str, Any]], selected: Collection[str]
) -> List[Dict[str, Any]]:
    """The rows whose status passes ``selected`` (all of them when it is empty)."""
    return [r for r in rows if status_selected(r.get("Validation"), selected)]


def consensus(validations: Iterable[Tuple[str, str]]) -> Tuple[str, str]:
    """Fold ``(status, inheritance)`` pairs into one ``(status, inheritance)``.

    Pass non-ignored validations only; an empty input gives ``("", "")``.

    - Different statuses give :data:`CONFLICTING`. ``in phase MNV`` agrees
      with ``present`` here, and wins over it in the result.
    - A present call takes an inheritance from its present validations:
      ``de novo`` first, then ``homozygous``, then the first one given.
    - ``absent`` stays ``absent``; anything else (``uncertain``, an SV's
      ``different``) is ``uncertain``.
    """
    pairs = [(s or "", i or "") for s, i in validations]
    if not pairs:
        return "", ""
    statuses = {s for s, _ in pairs}
    normalized = {"present" if s == IN_PHASE_MNV else s for s in statuses}
    if len(normalized) > 1:
        return CONFLICTING, ""
    if "present" in normalized:
        status = IN_PHASE_MNV if IN_PHASE_MNV in statuses else "present"
        inheritances = [i for s, i in pairs if s in _PRESENT and i]
        for preferred in ("de novo", "homozygous"):
            if preferred in inheritances:
                return status, preferred
        return status, inheritances[0] if inheritances else ""
    if "absent" in normalized:
        return "absent", ""
    return "uncertain", ""
