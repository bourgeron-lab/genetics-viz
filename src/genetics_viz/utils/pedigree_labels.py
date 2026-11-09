"""Display labels for the coded sex and phenotype values of a pedigree.

Pedigree files store sex and phenotype as bare codes, and the values reaching
the app are raw strings straight from the file - ``models.py`` reads pedigrees
with ``infer_schema_length=0`` and only normalises the parent columns. Some
files write them float-stringified (``1.0``, ``-9.0``, as the PMS cohort does),
so a plain comparison against ``"-9"`` misses them. Normalise through here
before grouping or labelling.
"""

from typing import Dict, List, Optional

from genetics_viz.models import is_dummy_id

#: Phenotype codes in display order: affected, unaffected, unknown.
PHENO_ORDER: List[str] = ["2", "1", "-9"]

PHENO_LABELS: Dict[str, str] = {
    "2": "2 (aff)",
    "1": "1 (unaff)",
    "-9": "-9 (unk)",
}

#: Sex codes in display order: male, female, unknown.
SEX_ORDER: List[str] = ["1", "2", "-9"]

SEX_LABELS: Dict[str, str] = {
    "1": "1 (male)",
    "2": "2 (female)",
    "-9": "-9 (unk)",
}

#: Plot colours per code. Phenotype matches the Statistics panel's palette.
PHENO_COLORS: Dict[str, str] = {"2": "#d62728", "1": "#2ca02c", "-9": "#7f7f7f"}
SEX_COLORS: Dict[str, str] = {"1": "#1f77b4", "2": "#e377c2", "-9": "#7f7f7f"}

_UNKNOWN = "-9"


def _strip_float(value: str) -> str:
    """Drop a trailing ``.0`` from a float-stringified code."""
    return value[:-2] if value.endswith(".0") else value


def normalize_pheno(pheno: Optional[str]) -> str:
    """Normalise a phenotype value to ``"1"``, ``"2"`` or ``"-9"``.

    Anything unrecognised - including ``None``, ``""`` and the ``"-"`` placeholder
    the cohort table uses - buckets under ``"-9"`` (unknown).
    """
    if not pheno:
        return _UNKNOWN
    pheno = _strip_float(pheno)
    return pheno if pheno in PHENO_LABELS else _UNKNOWN


def normalize_sex(sex: Optional[str]) -> str:
    """Normalise a sex value to ``"1"``, ``"2"`` or ``"-9"``."""
    if not sex:
        return _UNKNOWN
    sex = _strip_float(sex)
    return sex if sex in SEX_LABELS else _UNKNOWN


def pheno_label(pheno: Optional[str]) -> str:
    """Display label for a phenotype value, normalising it first."""
    return PHENO_LABELS[normalize_pheno(pheno)]


def sex_label(sex: Optional[str]) -> str:
    """Display label for a sex value, normalising it first."""
    return SEX_LABELS[normalize_sex(sex)]


# ---------------------------------------------------------------------------
# Dummy pedigree members (``models.is_dummy_id``)
# ---------------------------------------------------------------------------

#: The badge beside a dummy's ID, wherever members are listed.
DUMMY_LABEL = "dummy"

#: The badge's tooltip: what a dummy is, and why it has no data.
DUMMY_TOOLTIP = (
    "Dummy pedigree member: never sequenced, it only links relatives in the"
    " pedigree (shared parents of twins, grandparents of cousins). No data."
)

#: Muted grey, for a dummy's ID and badge.
DUMMY_COLOR = "#6b7280"
DUMMY_BACKGROUND = "#f3f4f6"


def dummy_badges(member_id: Optional[str]) -> List[Dict[str, str]]:
    """``badge_list`` cell badges for a member ID: one muted badge for a dummy.

    The badge is the ID itself, grey on grey, with :data:`DUMMY_TOOLTIP`.
    Empty for anyone else, so the cell shows the plain ID.
    """
    if not is_dummy_id(member_id):
        return []
    return [
        {
            "label": member_id,
            "color": DUMMY_BACKGROUND,
            "textColor": DUMMY_COLOR,
            "tooltip": DUMMY_TOOLTIP,
        }
    ]
