"""The family Diagnostics panel: one line per (kind, record, verdict).

Shared by the cohort family page and the standalone family page. Callers wrap
:func:`render_family_diagnostics` in their own ``@ui.refreshable``, so the
panel follows the member selection.
"""

from pathlib import Path
from typing import Dict, List

from nicegui import ui

from genetics_viz.components.diagnostic_loader import (
    format_diagnostic_contributors,
    load_family_diagnostics,
)
from genetics_viz.utils.cytobands import get_cytoband_range
from genetics_viz.utils.diagnostic_files import TYPE_COLORS
from genetics_viz.utils.diagnostic_labels import describe_str, describe_upd
from genetics_viz.utils.gene_scoring import (
    GeneScorer,
    build_gene_badges,
    gene_badges_html,
    get_gene_scorer,
    split_genes,
)

_DIAG_COLORS: Dict[str, str] = {
    "pathogenic": "red",
    "uncertain": "orange",
    "benign": "green",
    "conflicting": "amber",
}


def parse_sv_cytoband(variant: str) -> str:
    """Extract cytoband range from an SV variant key like 'chr1:1000-2000:dup'."""
    try:
        # Strip trailing :dup/:del if present
        parts = variant.split(":")
        if len(parts) >= 2:
            chrom = parts[0]
            range_part = parts[1]
            if "-" in range_part:
                start_str, end_str = range_part.split("-", 1)
                return get_cytoband_range(chrom, int(start_str), int(end_str))
    except (ValueError, IndexError):
        pass
    return ""


def render_gene_badges(gene_str: str, gene_scorer: GeneScorer) -> None:
    """Render the gene badges of a diagnostic entry, best tier first.

    Takes SV ``SYMBOL:type`` lists as well as plain symbols. Shows up to six
    genes, then a "+N genes" badge linking to all of them in GeneTrek.
    """
    symbols = [part.split(":")[0].strip() for part in split_genes(gene_str, ",")]
    badges = build_gene_badges([s for s in symbols if s and s != "-"], gene_scorer)
    ui.html(gene_badges_html(badges), sanitize=False)


def _kind_tag(kind: str) -> ui.label:
    return (
        ui.label(kind)
        .classes("text-xs text-white rounded px-1")
        .style(f"background-color: {TYPE_COLORS.get(kind, '#6b7280')}")
    )


def render_family_diagnostics(
    data_dir: Path, family_id: str, sample_ids: List[str]
) -> None:
    """The diagnostics recorded for ``sample_ids``, or a note that there are none."""
    entries = load_family_diagnostics(data_dir, family_id, sample_ids)
    if not entries:
        ui.label("No diagnostics recorded for selected members").classes(
            "text-gray-500 text-sm italic"
        )
        return

    gene_scorer = get_gene_scorer()
    for entry in entries:
        diag = entry.get("Diagnostic", "")
        source = entry.get("_source")
        variant = entry.get("Variant", "")
        gene_str = entry.get("Gene", "")
        with (
            ui.row()
            .classes("items-center gap-2 w-full px-2 py-1 border-b border-gray-100")
            .mark(f"diagnostic-entry diagnostic-{source}")
        ):
            ui.badge(
                diag.upper()[:3] if diag else "?",
                color=_DIAG_COLORS.get(diag, "grey"),
            ).classes("text-xs")

            if source == "str":
                # The gene first, then what was found: the key itself says little.
                if gene_str:
                    render_gene_badges(gene_str, gene_scorer)
                text, tooltip = describe_str(data_dir, family_id, entry)
                _kind_tag("STR").tooltip(tooltip)
                if text:
                    ui.label(text).classes("text-xs font-medium").tooltip(tooltip)
            elif source == "upd":
                text, tooltip = describe_upd(data_dir, family_id, entry)
                _kind_tag("UPD").tooltip(tooltip)
                ui.label(text).classes("text-xs font-medium").tooltip(tooltip)
            else:
                ui.label(variant).classes("text-xs font-mono")
                if source == "sv" and variant:
                    cytoband = parse_sv_cytoband(variant)
                    if cytoband:
                        ui.label(cytoband).classes("text-xs text-purple-600 italic")
                if gene_str:
                    render_gene_badges(gene_str, gene_scorer)

            ui.space()
            ui.label(entry.get("Sample", "")).classes("text-xs text-gray-600")
            # One row can merge several curators, so spell out who recorded
            # what on hover.
            user_label = ui.label(entry.get("User", "")).classes(
                "text-xs text-gray-400"
            )
            if len(entry.get("_entries", [])) > 1:
                user_label.tooltip(format_diagnostic_contributors(entry))
