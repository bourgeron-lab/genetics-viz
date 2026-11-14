"""The variant dialog's card of the denovo step's evidence for one candidate.

Shown when the dialog is opened from the De novo sub-tab with the step's
tables: the row carries the candidate's dnm row as ``_dnm``.
"""

from __future__ import annotations

from typing import Any, Dict

from nicegui import ui

from genetics_viz.utils import denovo

# ---------------------------------------------------------------------------


def _n(value: Any) -> str:
    if value is None or value == "":
        return "–"
    if isinstance(value, float):
        return f"{value:.3g}"
    return str(value)


def render_evidence_card(ev: Dict[str, Any]) -> None:
    """The denovo step's evidence for one candidate: reads, origin, flags."""
    flags = denovo.split_flags(ev.get("flags"))
    with ui.card().classes("w-full mb-4").mark("denovo-evidence"):
        with ui.column().classes("p-2 gap-2 w-full"):
            with ui.row().classes("items-center gap-2"):
                ui.label("De novo evidence").classes("text-lg font-semibold")
                strict = denovo.strict_badges(ev.get("strict"))
                for b in strict:
                    ui.badge(b["label"], color=b["color"]).tooltip(b["tooltip"])
                for b in denovo.flag_badges(flags):
                    ui.badge(b["label"], color=b["color"]).tooltip(b["tooltip"])
            reads = [
                ("", "Alt (HQ)", "Ref (HQ)", "Alt (all reads)", "Note"),
                (
                    f"Child {ev.get('child')}",
                    _n(ev.get("child_alt")),
                    _n(ev.get("child_ref")),
                    _n(ev.get("child_alt_raw")),
                    f"VAF {_n(ev.get('child_vaf'))}, binomial p {_n(ev.get('child_binom_p'))};"
                    f" alt fwd/rev {_n(ev.get('child_alt_fwd'))}/{_n(ev.get('child_alt_rev'))},"
                    f" ref fwd/rev {_n(ev.get('child_ref_fwd'))}/{_n(ev.get('child_ref_rev'))},"
                    f" strand p {_n(ev.get('strand_p'))}",
                ),
                (
                    f"Father {ev.get('father') or ''}",
                    _n(ev.get("father_alt")),
                    _n(ev.get("father_ref")),
                    _n(ev.get("father_alt_raw")),
                    "",
                ),
                (
                    f"Mother {ev.get('mother') or ''}",
                    _n(ev.get("mother_alt")),
                    _n(ev.get("mother_ref")),
                    _n(ev.get("mother_alt_raw")),
                    "",
                ),
            ]
            with (
                ui.element("div")
                .classes("grid gap-x-6 gap-y-0.5 text-sm")
                .style("grid-template-columns: repeat(5, max-content)")
            ):
                for i, line in enumerate(reads):
                    for cell in line:
                        ui.label(cell).classes(
                            "text-xs font-semibold text-gray-500" if i == 0 else ""
                        )
            with ui.row().classes("gap-6 flex-wrap text-sm"):
                ui.label(
                    "Siblings: "
                    + (
                        ev.get("sibs_with_alt")
                        or ev.get("sibs_called")
                        or "none with the call"
                    )
                ).tooltip("sample: high-quality alt reads / depth")
                origin = ev.get("origin")
                ui.label(
                    "Origin: "
                    + (
                        f"{origin} (evidence {_n(ev.get('origin_evidence'))})"
                        if origin
                        else denovo.origin_status_text(ev.get("origin_status")) or "–"
                    )
                )
                ui.label(
                    f"Cluster: {_n(ev.get('cluster_n'))} other candidate(s) within 20 bp,"
                    f" {_n(ev.get('nearby_calls'))} other call(s) within ±20 bp"
                )
                if ev.get("recurrence_families") not in (None, ""):
                    ui.label(f"Other families: {_n(ev.get('recurrence_families'))}")
                ui.label(
                    f"Wombat: GT {ev.get('wombat_child_gt') or '–'}, DP"
                    f" {_n(ev.get('wombat_child_dp'))}, GQ {_n(ev.get('wombat_child_gq'))},"
                    f" VAF {_n(ev.get('wombat_child_vaf'))}"
                ).classes("text-gray-600")
            if flags:
                with ui.column().classes("gap-0"):
                    for flag in sorted(flags, key=denovo.flag_rank):
                        info = denovo.flag_info(flag)
                        ui.label(f"{info['label']}: {info['tooltip']}").classes(
                            "text-xs text-gray-600"
                        )
