"""Pedigree drawings, through pypedigree (github.com/bourgeron-lab/pypedigree).

pypedigree lays the family out and writes the SVG; this module turns our
members into a pedigree text and picks the style. The output is
:func:`pypedigree.render_svg`: no script, every id under ``id_prefix``, CSS
scoped under ``.pypedigree``, so several drawings and the report's own SVGs
can share a page.

Black and white, standard notation: square, circle, diamond; affected
filled. Dummies (never sequenced, :func:`models.is_dummy_id`) have a dashed
outline (pypedigree's ``not_sampled`` column), as has the missing parent of a
member with one parent (pypedigree's phantom).

Things pypedigree 0.5.2 does not do, worked around here:
- **The text parser, not ``PedigreeBuilder``.** The builder leaves a child
  with one declared parent unconnected; the parser gives it a phantom
  second parent (bourgeron-lab/pypedigree#40).
- **A parent named but absent from the list is written as a dashed row**,
  rather than left to the parser, which would draw the same phantom but warn,
  and ``warnings.catch_warnings`` is not thread-safe (the report is built in
  ``asyncio.to_thread``).
- **The spacing follows the longest label.** The layout spaces symbols, not
  labels, so ``DUMMY_<barcode>_MF`` beside its partner's overlapped
  (bourgeron-lab/pypedigree#41).

Pure and fast (about 1 ms for a trio); ``import pypedigree`` does not load
matplotlib.
"""

import logging
import re
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import pypedigree

logger = logging.getLogger(__name__)

#: Label size in points. pypedigree's default (8) is small beside a 0.8 symbol.
FONT_SIZE = 10.0
#: One point of font in layout units (pypedigree draws 8 pt at 0.16).
_UNITS_PER_POINT = 0.02
#: Average width of a label character, in ems (Helvetica, mostly capitals and
#: digits).
_CHAR_EM = 0.66
#: Clear space between two neighbouring labels, in layout units.
_LABEL_GAP = 0.3

#: Our normalised pedigree codes (:mod:`pedigree_labels`) -> PED codes.
_SEX = {"1": "1", "2": "2"}
_PHENOTYPE = {"1": "1", "2": "2"}

_HEADER = (
    "FID",
    "IID",
    "Father",
    "Mother",
    "Sex",
    "Phenotype",
    "not_sampled",
    "a_Note",
)
_VIEWBOX = re.compile(r'viewBox="([-\d.]+) ([-\d.]+) ([\d.]+) ([\d.]+)"')
_SPACE = re.compile(r"\s+")


@dataclass
class PedigreePerson:
    """One individual of the drawing; codes as :mod:`pedigree_labels` writes them."""

    iid: str
    father: str = ""
    mother: str = ""
    #: ``"1"`` male, ``"2"`` female, anything else unknown.
    sex: str = "-9"
    #: ``"2"`` affected, ``"1"`` unaffected, anything else unknown.
    phenotype: str = "-9"
    #: Never sequenced: dashed outline, status unknown.
    dummy: bool = False
    #: A line written under the ID (for instance an atypical karyotype).
    note: str = ""


@dataclass
class PedigreeDrawing:
    svg: str
    #: The ``viewBox`` size, in pypedigree's layout units (a symbol is 0.8).
    width: float
    height: float
    #: Someone is drawn dashed: a dummy, or a parent who is not sequenced.
    has_unsampled: bool = False
    #: Someone has a note under their ID.
    has_notes: bool = False


def _cell(value: str) -> str:
    """One PED cell: whitespace would split it, an empty one means missing."""
    return _SPACE.sub("_", value.strip())


def _absent_parents(people: Sequence[PedigreePerson]) -> List[Tuple[str, str]]:
    """``(id, sex code)`` of every parent named but not in ``people``, in order."""
    known = {p.iid for p in people}
    absent: List[Tuple[str, str]] = []
    for p in people:
        for parent, sex in ((p.father, "1"), (p.mother, "2")):
            if parent and parent not in known:
                known.add(parent)
                absent.append((parent, sex))
    return absent


def pedigree_text(family_id: str, people: Sequence[PedigreePerson]) -> str:
    """The PED text pypedigree parses: our members, then the absent parents."""
    fid = _cell(family_id) or "family"
    lines = ["\t".join(_HEADER)]
    for p in people:
        lines.append(
            "\t".join(
                [
                    fid,
                    _cell(p.iid),
                    _cell(p.father) or "0",
                    _cell(p.mother) or "0",
                    _SEX.get(p.sex, "0"),
                    "0" if p.dummy else _PHENOTYPE.get(p.phenotype, "0"),
                    "1" if p.dummy else "0",
                    _cell(p.note),
                ]
            )
        )
    for parent, sex in _absent_parents(people):
        lines.append("\t".join([fid, _cell(parent), "0", "0", sex, "0", "1", ""]))
    return "\n".join(lines) + "\n"


def _style(people: Sequence[PedigreePerson]) -> "pypedigree.PedigreeStyle":
    """pypedigree's style, spaced so the longest label clears its neighbours."""
    style = pypedigree.PedigreeStyle()
    style.font_size = FONT_SIZE
    longest = max(len(_cell(p.iid)) for p in people)
    for parent, _ in _absent_parents(people):
        longest = max(longest, len(_cell(parent)))
    label = longest * FONT_SIZE * _UNITS_PER_POINT * _CHAR_EM
    spacing = label + _LABEL_GAP
    style.couple_gap = max(style.couple_gap, spacing)
    style.sibling_gap = max(style.sibling_gap, spacing)
    style.min_node_distance = max(style.min_node_distance, spacing)
    # A label wider than its symbol overhangs it on both sides.
    style.margin = max(style.margin, (label - style.shape_size) / 2 + 0.15)
    return style


def pedigree_svg(
    family_id: str,
    people: Sequence[PedigreePerson],
    *,
    id_prefix: str = "ped-",
) -> Optional[PedigreeDrawing]:
    """The family drawn as a static, responsive SVG, ``None`` if it cannot be.

    A pedigree pypedigree refuses is logged and gives ``None``: the report
    keeps its table, as for any other missing plot.
    """
    if not people:
        return None
    try:
        pedigree = pypedigree.parse_pedigree_from_string(
            pedigree_text(family_id, people)
        )
        svg = pypedigree.render_svg(
            pedigree,
            style=_style(people),
            show_legend=False,
            id_prefix=id_prefix,
            responsive=True,
        )
    except Exception:
        logger.warning("Failed to draw the pedigree of %s", family_id, exc_info=True)
        return None
    match = _VIEWBOX.search(svg)
    width, height = (float(match.group(3)), float(match.group(4))) if match else (0, 0)
    return PedigreeDrawing(
        svg=svg,
        width=width,
        height=height,
        has_unsampled=bool(_absent_parents(people))
        or any(p.dummy or bool(p.father) != bool(p.mother) for p in people),
        has_notes=any(p.note for p in people),
    )
