"""SV dialog for a Vedette call: its details card and :class:`SVSource`.

The card shows the family site (type, size, impact, genes, cohort frequency,
every member's genotype) and the carrier's own call before merging, when it
has one: its FILTER, evidence source, breakpoint precision and read counts.
Validations go to ``validations/vedette.tsv`` with the site id, the caller
version and the model version (``VEDETTE_VALIDATION_HEADER``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from nicegui import ui

from genetics_viz.components.sv_dialog import SVSource
from genetics_viz.components.validation_loader import VEDETTE_VALIDATION_HEADER
from genetics_viz.utils.gene_scoring import (
    build_gene_badges,
    gene_badges_html,
    get_gene_scorer,
    tier_legend_html,
)
from genetics_viz.utils.genetrek import search_url
from genetics_viz.utils.vedette import (
    OwnCall,
    assume_precise,
    find_own_call,
    sample_vedette_vcf,
)

_FIELD_LABELS = {
    "GT": "GT",
    "GQ": "GQ",
    "SQ": "SQ (calibrated quality)",
    "CN": "CN",
    "CNQ": "CN posterior",
    "RD": "Copy ratio",
    "PE": "Discordant pairs",
    "SR": "Split / clipped reads",
    "CG": "Reads with the CIGAR deletion",
}


def _badge(text: str, color: str, text_color: str = "white") -> str:
    return (
        f'<span class="q-badge" style="background-color: {color}; color: '
        f"{text_color}; padding: 2px 6px; border-radius: 4px; font-size: 12px;"
        f'">{text}</span>'
    )


def _field(label: str, value: Any) -> None:
    with ui.column().classes("gap-0"):
        ui.label(label).classes("text-xs font-semibold text-gray-500")
        ui.label("-" if value in (None, "") else str(value)).classes(
            "text-base text-gray-900 font-medium"
        )


def _fmt_af(value: Optional[float]) -> str:
    return "-" if value is None else f"{value:.4g}"


def render_vedette_details(row: Dict[str, Any], own: Optional[OwnCall]) -> None:
    """The "SV Details" card of a Vedette carrier row."""
    ui.label("SV Details (Vedette)").classes("text-xl font-semibold mb-2")
    scorer = get_gene_scorer()
    with ui.card().classes("w-full mb-4"):
        with ui.column().classes("p-4 gap-4 w-full"):
            with ui.row().classes("gap-4 items-center flex-wrap"):
                ui.html(tier_legend_html(scorer), sanitize=False)
                ui.link(
                    "Region in GeneTrek",
                    search_url([row.get("chr:start-end", "")]),
                    new_tab=True,
                ).classes("text-xs")

            ui.separator()
            with ui.row().classes("gap-6 flex-wrap items-center"):
                with ui.column().classes("gap-0"):
                    ui.label("type").classes("text-xs font-semibold text-gray-500")
                    ui.html(
                        _badge(row.get("svtype", ""), row.get("_type_color", "#999")),
                        sanitize=False,
                    )
                _field("size (bp)", f"{row.get('size', 0):,}")
                _field("carrier", row.get("sample"))
                _field(
                    "GT · CN · GQ",
                    f"{row.get('GT')} · {row.get('CN')} · {row.get('GQ')}",
                )
                with ui.column().classes("gap-0"):
                    ui.label("inheritance").classes(
                        "text-xs font-semibold text-gray-500"
                    )
                    if row.get("inheritance"):
                        ui.html(
                            _badge(row["inheritance"], row["_inheritance_color"]),
                            sanitize=False,
                        )
                    else:
                        ui.label("-").classes("text-base")
                with ui.column().classes("gap-0"):
                    ui.label("impact").classes("text-xs font-semibold text-gray-500")
                    with ui.row().classes("gap-2 items-center"):
                        if row.get("impact_level"):
                            ui.html(
                                _badge(row["impact_level"], row["_impact_color"]),
                                sanitize=False,
                            )
                        score = row.get("impact_score")
                        if score is not None:
                            ui.label(f"{score:+.2f}").classes("text-base font-medium")
                _field("cohort AF", _fmt_af(row.get("cohort_af")))
                _field("carrier fraction", _fmt_af(row.get("cohort_carrier_frac")))
                _field("carriers", row.get("cohort_ncarrier"))
                _field("site", row.get("site_id"))

            if row.get("impact_reason"):
                ui.label(row["impact_reason"]).classes("text-sm text-gray-700")
            genes = row.get("_genes") or []
            if genes:
                with ui.column().classes("gap-1"):
                    ui.label("genes").classes("text-xs font-semibold text-gray-500")
                    ui.html(
                        gene_badges_html(build_gene_badges(genes, scorer, cap=None)),
                        sanitize=False,
                    )
                    if "(+" in (row.get("genes") or ""):
                        ui.label(
                            "Vedette truncated the gene list: the region link "
                            "above lists every gene."
                        ).classes("text-xs text-gray-500 italic")

            with ui.column().classes("gap-1"):
                ui.label("family (GT·CN)").classes(
                    "text-xs font-semibold text-gray-500"
                )
                with ui.row().classes("gap-2 flex-wrap"):
                    for b in row.get("family_badges", []):
                        ui.html(
                            _badge(b["label"], b["color"], b["textColor"]),
                            sanitize=False,
                        ).tooltip(b["tooltip"])

            ui.separator()
            _render_own_call(row, own)


def _render_own_call(row: Dict[str, Any], own: Optional[OwnCall]) -> None:
    sample = row.get("sample", "")
    if own is None:
        ui.label(
            f"{sample} has no own call of this type here: re-genotyped at the "
            "cohort site."
        ).classes("text-sm text-gray-600 italic").mark("vedette-own-call")
        return
    with ui.column().classes("gap-2").mark("vedette-own-call"):
        with ui.row().classes("items-center gap-3 flex-wrap"):
            ui.label(f"{sample}'s own call").classes("text-sm font-semibold")
            ui.label(f"{own.call_id} · {own.start:,}-{own.end:,}").classes(
                "text-sm font-mono"
            )
            filters = ", ".join(own.filters) or "PASS"
            ui.html(
                _badge(filters, "#16a34a" if own.passed else "#9ca3af"),
                sanitize=False,
            )
            if own.source:
                ui.html(_badge(own.source, "#475569"), sanitize=False)
            ui.html(
                _badge(
                    "precise" if own.precise else "imprecise",
                    "#0f766e" if own.precise else "#a16207",
                ),
                sanitize=False,
            )
            if own.model:
                ui.label(f"model {own.model}").classes("text-xs text-gray-500")
        with ui.row().classes("gap-6 flex-wrap"):
            for key, label in _FIELD_LABELS.items():
                if key in own.fields:
                    value = own.fields[key]
                    if isinstance(value, float):
                        value = f"{value:.3g}"
                    _field(label, value)
            if own.cipos:
                _field("CIPOS", f"{own.cipos[0]},{own.cipos[1]}")
            if own.ciend:
                _field("CIEND", f"{own.ciend[0]},{own.ciend[1]}")


def vedette_source(
    data_dir: Path, family_id: str, row: Dict[str, Any], caller_version: str
) -> SVSource:
    """The dialog's source for a carrier row: own call read here (one tabix query)."""
    own = find_own_call(
        sample_vedette_vcf(data_dir, row["sample"]),
        row["contig"],
        int(row["start"]),
        int(row["end"]),
        row["svtype"],
    )
    precise = own.precise if own is not None else assume_precise(row.get("size"))
    return SVSource(
        validation_kind="Vedette",
        header=list(VEDETTE_VALIDATION_HEADER),
        extra_fields={
            "SiteId": row.get("site_id", ""),
            "CallerVersion": caller_version,
            "ModelVersion": own.model if own is not None else "",
        },
        precise=precise,
        render_details=lambda: render_vedette_details(row, own),
        label="Vedette",
    )
