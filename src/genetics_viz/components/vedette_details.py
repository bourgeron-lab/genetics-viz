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
    ExportInfo,
    OwnCall,
    resolve_own_call,
    row_is_precise,
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


_OWN_ORIGIN = {
    "parquet": "from Vedette's export",
    "tsv": "from Vedette's export (code, SQ and source only)",
    "vcf": "matched in the sample's VCF",
}


def _fmt_af(value: Optional[float]) -> str:
    return "-" if value is None else f"{value:.4g}"


def render_vedette_details(row: Dict[str, Any], own: Optional[OwnCall]) -> None:
    """The "SV Details" card of a Vedette carrier row."""
    ui.label("SV Details (Vedette)").classes("text-xl font-semibold mb-2")
    scorer = get_gene_scorer()
    with ui.card().classes("w-full mb-4"):
        with ui.column().classes("p-4 gap-4 w-full"):
            with ui.row().classes("gap-4 items-center flex-wrap"):
                ui.html(tier_legend_html(scorer, exonic_border=True), sanitize=False)
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
                for label, code_field, color_field in (
                    ("kind", "kind", "_kind_color"),
                    ("own call", "own_call", "_own_call_color"),
                    ("CN support", "support", "_support_color"),
                ):
                    if row.get(code_field):
                        with ui.column().classes("gap-0"):
                            ui.label(label).classes(
                                "text-xs font-semibold text-gray-500"
                            )
                            ui.html(
                                _badge(row[code_field], row[color_field]),
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
                        dnq = row.get("DNQ")
                        ui.html(
                            _badge(
                                row["inheritance"]
                                + (f" · DNQ {dnq}" if dnq is not None else ""),
                                row["_inheritance_color"],
                            ),
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
                _field(
                    "cohort carrier fraction", _fmt_af(row.get("cohort_carrier_frac"))
                )
                if row.get("ref_carrier_frac") is not None:
                    _field(
                        "reference carrier fraction", _fmt_af(row["ref_carrier_frac"])
                    )
                if row.get("cohort_af") is not None:  # exports before v1.4.0
                    _field("cohort AF", _fmt_af(row["cohort_af"]))
                if row.get("sqmax") is not None:
                    _field("SQ max", f"{row['sqmax']:g}")
                if row.get("mosaic_cn") is not None:
                    _field("mosaic CN", f"{row['mosaic_cn']:g}")
                for label, key in (
                    ("LCR overlap", "lcr_ovl"),
                    ("segdup overlap", "segdup_ovl"),
                ):
                    if row.get(key):
                        _field(label, f"{row[key]:.0%}")
                flt = row.get("filter") or ""
                if flt and flt != "PASS":
                    with ui.column().classes("gap-0"):
                        ui.label("filter").classes(
                            "text-xs font-semibold text-gray-500"
                        )
                        ui.html(_badge(flt, "#9ca3af"), sanitize=False)
                else:
                    _field("filter", flt)
                for label, key, fmt in (
                    ("cohort AF", "af", ".4g"),
                    ("samples merged", "nsamples", "d"),
                    ("precise fraction", "precise_frac", ".0%"),
                    ("SQ mean", "sqmean", ".3g"),
                    ("HWE p", "hwe_p", ".2g"),
                    ("CN discordance", "cn_disc", ".0%"),
                    ("reference AF", "ref_af", ".4g"),
                ):
                    if row.get(key) is not None:
                        _field(label, format(row[key], fmt))
                if row.get("rescue"):
                    _field("rescue", row["rescue"])
                ids = row.get("site_ids") or []
                _field(
                    "cohort sites" if len(ids) > 1 else "site",
                    ", ".join(ids) or row.get("site_id"),
                )
                if row.get("site_key"):
                    _field("site key", row["site_key"])

            if row.get("qc_alerts"):
                ui.label(f"Vedette alerts: {row['qc_alerts']}").classes(
                    "text-sm text-amber-800"
                )

            if row.get("impact_reason"):
                ui.label(row["impact_reason"]).classes("text-sm text-gray-700")
            genes = row.get("_genes") or []
            if genes:
                with ui.column().classes("gap-1"):
                    ui.label("genes").classes("text-xs font-semibold text-gray-500")
                    ui.html(
                        gene_badges_html(
                            build_gene_badges(
                                genes,
                                scorer,
                                cap=None,
                                bordered=row.get("_exonic") or (),
                            )
                        ),
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
            f"{sample} has no own call of this type here: carried through the "
            "cohort's re-genotyping only."
        ).classes("text-sm text-gray-600 italic").mark("vedette-own-call")
        return
    with ui.column().classes("gap-2").mark("vedette-own-call"):
        with ui.row().classes("items-center gap-3 flex-wrap"):
            ui.label(f"{sample}'s own call").classes("text-sm font-semibold")
            where = " · ".join(
                part
                for part in (
                    own.call_id,
                    f"{own.start:,}-{own.end:,}"
                    if own.start is not None and own.end is not None
                    else "",
                )
                if part
            )
            if where:
                ui.label(where).classes("text-sm font-mono")
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
            ui.label(_OWN_ORIGIN.get(own.origin, "")).classes(
                "text-xs text-gray-400 italic"
            )
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
    data_dir: Path, family_id: str, row: Dict[str, Any], info: ExportInfo
) -> SVSource:
    """The dialog's source for a carrier row.

    The own call comes from the export when it carries it (Vedette v1.9.0),
    else from one tabix query of the sample's VCF. Validations are keyed by
    the row's ``_variant_key`` (the cohort site's ``site_key`` when there is
    one) and record the row id, the caller version and the model version: the
    export's (JSON), else the own call's ``MODEL_VER``.
    """
    own = resolve_own_call(data_dir, row)
    precise = row_is_precise(row, own)
    return SVSource(
        validation_kind="Vedette",
        header=list(VEDETTE_VALIDATION_HEADER),
        extra_fields={
            "SiteId": row.get("site_id", ""),
            "CallerVersion": info.caller,
            "ModelVersion": info.model_version
            or (own.model if own is not None else ""),
        },
        precise=precise,
        render_details=lambda: render_vedette_details(row, own),
        label="Vedette",
        # Vedette's start is the first affected base: curated starts match it.
        start_offset=1,
        variant_key=row.get("_variant_key"),
    )
