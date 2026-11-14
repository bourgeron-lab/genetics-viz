"""The standing caveat on a cell-line sample's results (ghfc-ngs ``genomeqc``).

A genome from an EBV-transformed cell line has sound inherited genotypes, but
its de novo calls, isodisomies and CNVs may come from the culture. Wherever
one of those is reviewed (UPD and SV dialogs and tabs, the diagnostic
dialog) a short note says so, acknowledged flag or not. The dialogs open from
pages without a cohort (Search, the validation pages), so the note finds the
cohorts holding the family itself.
"""

import logging
from typing import Any, Dict, Iterable, List, Optional

from nicegui import ui

from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.genome_qc import caveat_text, cell_lines
from genetics_viz.utils.sample_qc import FlagRow

logger = logging.getLogger(__name__)


def family_cohorts(family_id: str, cohort: Any = None) -> List[Any]:
    """``[cohort]``, or every cohort of the user's data directory holding the family."""
    if cohort is not None:
        return [cohort]
    try:
        store = get_data_store()
    except Exception:
        return []
    return [
        c for c in store.cohorts.values() if family_id in getattr(c, "families", {})
    ]


def family_cell_lines(family_id: str, cohort: Any = None) -> Dict[str, FlagRow]:
    """Sample -> ``CELL_LINE`` flag of the family's members; ``{}`` on any failure."""
    try:
        return cell_lines(family_cohorts(family_id, cohort), family_id)
    except Exception:
        logger.warning("Failed to read the genome QC of %s", family_id, exc_info=True)
        return {}


def render_cell_line_note(
    family_id: str,
    samples: Iterable[str],
    context: str,
    cohort: Any = None,
    known: Optional[Dict[str, FlagRow]] = None,
) -> bool:
    """The caveat for the cell lines among ``samples``; whether one was drawn.

    ``context`` is a key of ``caveats`` in ``genome_qc.yaml`` (``upd``,
    ``cnv``, ``diagnostic``). ``known`` skips the lookup.
    """
    found = known if known is not None else family_cell_lines(family_id, cohort)
    text = caveat_text(context, [s for s in samples if s in found])
    if not text:
        return False
    with (
        ui.row()
        .classes(
            "w-full items-start gap-2 p-2 rounded border border-red-300 bg-red-50"
            " no-wrap"
        )
        .mark(f"cell-line-note-{context}")
    ):
        ui.icon("biotech", color="red-8").classes("mt-0.5")
        ui.label(text).classes("text-sm text-red-900")
    return True
