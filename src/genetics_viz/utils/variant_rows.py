"""Wombat variant rows: their badges, and the filters the tables apply to them.

The Wombat tab and the Search page show the same rows with the same filters
(consequences, LCR, gnomAD, MNV); both go through here. The functions are
pure: the pages own the filter state and the table.

Aggregated rows (``wombat.aggregate_by_variant``) join a column's values with
``,``, and VEP joins several terms of one transcript with ``&``, so a cell is
split on both.
"""

from typing import Any, Collection, Dict, Iterable, List, Set

from genetics_viz.utils.clinvar import (
    CLINVAR_COLORS,
    format_clinvar_display,
    get_clinvar_color,
)
from genetics_viz.utils.score_colors import get_score_color
from genetics_viz.utils.vep import (
    VEP_CONSEQUENCES,
    format_consequence_display,
    get_consequence_color,
)

_CLINVAR_KNOWN = {k.lower() for k in CLINVAR_COLORS}

#: An ``mnv_proba`` above this marks the call as part of an MNV.
MNV_PROBA_THRESHOLD = 0.5


def split_terms(value: Any) -> List[str]:
    """The terms of a VEP cell: split on ``,`` and ``&``, blanks and ``.`` dropped."""
    if not value:
        return []
    terms = []
    for part in str(value).split(","):
        for term in part.split("&"):
            term = term.strip()
            if term and term != ".":
                terms.append(term)
    return terms


def _badges(terms: Iterable[str], label_of, color_of) -> List[Dict[str, str]]:
    """One badge per distinct (label, colour), in first-seen order."""
    badges: List[Dict[str, str]] = []
    seen = set()
    for term in terms:
        key = (label_of(term), color_of(term))
        if key not in seen:
            seen.add(key)
            badges.append({"label": key[0], "color": key[1]})
    return badges


def consequence_badges(value: Any, unknown: Set[str]) -> List[Dict[str, str]]:
    """Badges for a ``VEP_Consequence`` cell; terms missing from the YAML go to ``unknown``."""
    terms = split_terms(value)
    unknown.update(t for t in terms if t not in VEP_CONSEQUENCES)
    return _badges(terms, format_consequence_display, get_consequence_color)


def clinvar_badges(value: Any, unknown: Set[str]) -> List[Dict[str, str]]:
    """Badges for a ClinVar significance cell; unknown terms go to ``unknown``."""
    terms = split_terms(value)
    unknown.update(t for t in terms if t.lower() not in _CLINVAR_KNOWN)
    return _badges(terms, format_clinvar_display, get_clinvar_color)


def add_vep_badges(
    row: Dict[str, Any], unknown_consequences: Set[str], unknown_clinvar: Set[str]
) -> None:
    """Set ``ConsequenceBadges`` and ``ClinVarBadges`` from the VEP columns."""
    row["ConsequenceBadges"] = consequence_badges(
        row.get("VEP_Consequence", ""), unknown_consequences
    )
    row["ClinVarBadges"] = clinvar_badges(row.get("VEP_CLIN_SIG", ""), unknown_clinvar)


def add_score_badges(row: Dict[str, Any]) -> None:
    """Add ``<column>_badge`` for every column with a continuous score colour scale."""
    # list(): the loop adds keys to the row it walks.
    for col_name, value_str in list(row.items()):
        if not value_str or value_str == ".":
            continue
        try:
            value = float(value_str)
        except (ValueError, TypeError):
            continue
        badge_info = get_score_color(col_name, value)
        if badge_info:
            row[f"{col_name}_badge"] = {
                "label": f"{value:.3f}",
                "color": badge_info["color"],
                "tooltip": f"{col_name}: {value:.3f} ({badge_info['label']})",
            }


def in_lcr(row: Dict[str, Any]) -> bool:
    """The call lies in a low-complexity region."""
    return bool(row.get("LCR")) and "true" in str(row.get("LCR", "")).lower()


def gnomad_filtered(row: Dict[str, Any]) -> bool:
    """gnomAD genomes flag the site (``genomes_filters`` is set)."""
    return bool(row.get("genomes_filters"))


def in_gnomad_wgs(row: Dict[str, Any]) -> bool:
    """The site has a gnomAD genomes filtering allele frequency."""
    return bool(row.get("fafmax_faf95_max_genomes"))


def likely_mnv(row: Dict[str, Any]) -> bool:
    """Any ``mnv_proba`` value of the row is above :data:`MNV_PROBA_THRESHOLD`."""
    for part in str(row.get("mnv_proba") or "").split(","):
        try:
            if float(part) > MNV_PROBA_THRESHOLD:
                return True
        except ValueError:
            continue
    return False


def consequences_selected(row: Dict[str, Any], selected: Collection[str]) -> bool:
    """Whether the row has a consequence in ``selected``.

    An empty selection, or every known consequence, filters nothing. A row
    with no consequence passes only then.
    """
    if not selected or set(selected) >= set(VEP_CONSEQUENCES):
        return True
    return any(t in selected for t in split_terms(row.get("VEP_Consequence", "")))
