"""The De novo sub-tab's own parts: the flag filter, the per-child summary, the
banner of curated calls the filters hide Data and wording come from
``utils/denovo.py``; the variant dialog's evidence card is
``components/denovo_evidence.py``.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Iterable, List, Optional

from nicegui import ui

from genetics_viz.utils import denovo
from genetics_viz.utils.denovo import DenovoData, DenovoFilter

#: Link keys of the de novo filter: mode, hidden flags, required flags.
DENOVO_LINK_KEYS = ("wdm", "wfh", "wfo")

_MODES = {"strict": "Strict", "all": "All candidates"}


def filter_from_link(link: Any) -> DenovoFilter:
    f = DenovoFilter()
    if link is None:
        return f
    mode = link.get("wdm")
    if mode in _MODES:
        f.mode["value"] = mode
    for key, state in (("wfh", f.hide), ("wfo", f.only)):
        values = link.get_list(key)
        if values:
            state["value"] = [v for v in values if v]
    return f


def filter_to_link(f: Optional[DenovoFilter], link: Any) -> None:
    if f is None:
        link.drop(*DENOVO_LINK_KEYS)
        return
    link.set("wdm", f.mode["value"], default="strict")
    link.set("wfh", f.hide["value"], default=[])
    link.set("wfo", f.only["value"], default=[])


def _badge_html(flag: str) -> str:
    info = denovo.flag_info(flag)
    return (
        f'<span class="dt-badge" style="background-color:{info["color"]};'
        f'color:white">{info["label"]}</span>'
    )


def render_denovo_filter(
    f: DenovoFilter, present: Iterable[str], on_change: Callable[[], Any]
) -> None:
    """Strict / All, and a menu setting each flag present to –, Hide or Only."""
    flags = sorted(set(present), key=denovo.flag_rank)

    @ui.refreshable
    def summary() -> None:
        parts = []
        if f.hide["value"]:
            parts.append(f"hiding {len(f.hide['value'])}")
        if f.only["value"]:
            parts.append(f"only {len(f.only['value'])}")
        ui.label("Flags" + (f" ({', '.join(parts)})" if parts else "")).classes(
            "normal-case"
        )

    def set_mode(e: Any) -> None:
        f.mode["value"] = e.value
        on_change()

    def make_flag_handler(flag: str) -> Callable[[Any], None]:
        def handler(e: Any) -> None:
            for state in (f.hide, f.only):
                if flag in state["value"]:
                    state["value"] = [x for x in state["value"] if x != flag]
            if e.value in ("hide", "only"):
                target = f.hide if e.value == "hide" else f.only
                target["value"] = target["value"] + [flag]
            if (
                e.value == "only"
                and f.mode["value"] == "strict"
                and flag not in denovo.STRICT_IGNORED
                and flag != denovo.NOT_ASSESSED
            ):
                f.mode["value"] = "all"
                mode_toggle.value = "all"
                ui.notify("Showing all candidates: Strict hides that flag")
            summary.refresh()
            on_change()

        return handler

    with ui.row().classes("items-center gap-3 mt-2 flex-wrap").mark("denovo-filter"):
        ui.label("De novo candidates:").classes("font-semibold text-blue-700")
        mode_toggle = (
            ui.toggle(_MODES, value=f.mode["value"], on_change=set_mode)
            .props("dense no-caps toggle-color=blue")
            .mark("denovo-mode")
        )
        mode_toggle.tooltip(
            "Strict: no flag but Sibling shared (and candidates the step did not"
            " assess). All: every candidate of the Wombat de novo table."
        )
        with ui.button(icon="flag").props("outline color=blue size=sm no-caps"):
            summary()
            with ui.menu().props("auto-close=false"), ui.column().classes("p-3 gap-1"):
                ui.label("Hide or keep only candidates with a flag").classes(
                    "text-xs text-gray-500"
                )
                for flag in flags:
                    value = (
                        "hide"
                        if flag in f.hide["value"]
                        else "only"
                        if flag in f.only["value"]
                        else ""
                    )
                    with ui.row().classes("items-center gap-2 no-wrap w-full"):
                        ui.html(_badge_html(flag), sanitize=False).classes(
                            "min-w-[9rem]"
                        ).tooltip(f"{flag}: {denovo.flag_info(flag)['tooltip']}")
                        ui.toggle(
                            {"": "–", "hide": "Hide", "only": "Only"},
                            value=value,
                            on_change=make_flag_handler(flag),
                        ).props("dense no-caps size=sm").mark(f"denovo-flag-{flag}")


def render_child_summary(data: DenovoData, members: List[str]) -> None:
    """One line per child: candidates, strict, rate, QC, origins."""
    rows = [q for q in data.qc_rows() if q.get("child") in members] or data.qc_rows()
    if not rows:
        return
    with (
        ui.element("div")
        .classes("w-full mt-2 rounded border border-gray-200 bg-gray-50 p-2")
        .mark("denovo-children")
    ):
        with ui.row().classes("items-center gap-2"):
            ui.label(f"denovo step · {data.label}").classes(
                "text-sm font-semibold text-gray-700"
            )
            masks = data.meta.get("masks")
            if masks:
                ui.label(f"masks {masks}").classes("text-xs text-gray-500")
            if data.level != "cohort":
                ui.label("(family tables: no cohort recurrence or QC verdict)").classes(
                    "text-xs text-gray-500"
                )
        header = [
            "Child", "Candidates", "Strict (SNV / indel)", "Strict per Gb",
            "QC", "Paternal / maternal", "Sibling shared",
        ]  # fmt: skip
        with (
            ui.element("div")
            .classes("grid gap-x-4 gap-y-0.5 text-sm mt-1")
            .style("grid-template-columns: repeat(7, max-content)")
        ):
            for h in header:
                ui.label(h).classes("text-xs font-semibold text-gray-500")
            for q in rows:
                ui.label(str(q.get("child")))
                ui.label(str(q.get("n_candidates") or 0))
                ui.label(
                    f"{q.get('n_strict') or 0}"
                    f" ({q.get('n_strict_snv') or 0} / {q.get('n_strict_indel') or 0})"
                ).tooltip("Autosomes only, as the per-child QC counts them")
                rate = q.get("rate_per_gb")
                ui.label(f"{rate:.1f}" if isinstance(rate, (int, float)) else "–")
                if q.get("qc"):
                    info = denovo.qc_info(q.get("qc"))
                    low, high = q.get("qc_low"), q.get("qc_high")
                    tip = info["tooltip"]
                    if isinstance(low, (int, float)) and isinstance(high, (int, float)):
                        tip += (
                            f" Range {low:.1f}–{high:.1f} per Gb, from"
                            f" {denovo.qc_reference_text(q.get('qc_reference'))}."
                        )
                    if q.get("sampleqc_flags"):
                        tip += f" Sample QC: {q.get('sampleqc_flags')}."
                    ui.badge(info["label"], color=info["color"]).tooltip(tip).classes(
                        "w-fit"
                    ).mark(f"denovo-qc-{q.get('child')}")
                else:
                    ui.label("–").tooltip("The QC verdict is in the cohort table")
                ui.label(f"{q.get('n_paternal') or 0} / {q.get('n_maternal') or 0}")
                ui.label(str(q.get("n_sib_shared") or 0))


def _row_text(r: Dict[str, Any]) -> str:
    flags = denovo.split_flags(r.get("dnm_flags") or r.get("flags"))
    labels = ", ".join(denovo.flag_info(f)["label"] for f in flags)
    what = r.get("Validation") or r.get("validation") or ""
    diag = r.get("Diagnostic") or r.get("diagnostic")
    return (
        f"{r.get('Variant')} in {r.get('sample')}"
        + (f" — {what}" if what else "")
        + (f", diagnostic {diag}" if diag else "")
        + (f" (flags: {labels})" if labels else "")
    )


def render_curation_banner(
    hidden: List[Dict[str, Any]],
    missing: List[Dict[str, Any]],
    on_show_all: Optional[Callable[[], None]],
    limit: int = 6,
) -> None:
    """Curated de novo calls the filters hide, and those no longer candidates."""
    if not hidden and not missing:
        return
    with (
        ui.row()
        .classes(
            "w-full items-start gap-2 p-3 mt-2 rounded border border-amber-400"
            " bg-amber-50 no-wrap"
        )
        .mark("denovo-curation-banner")
    ):
        ui.icon("visibility_off", color="amber-9").classes("mt-0.5")
        with ui.column().classes("gap-0 min-w-0 flex-1"):
            if hidden:
                ui.label(
                    f"{len(hidden)} curated de novo call(s) hidden by the current"
                    " filters (validated present without an inherited mode, or"
                    " with a diagnostic):"
                ).classes("font-semibold text-amber-900")
                for r in hidden[:limit]:
                    ui.label(f"• {_row_text(r)}").classes("text-sm text-amber-900")
                if len(hidden) > limit:
                    ui.label(f"… and {len(hidden) - limit} more").classes(
                        "text-sm text-amber-900"
                    )
            if missing:
                ui.label(
                    f"{len(missing)} validated de novo call(s) of this family are not"
                    " candidates in the Wombat de novo table:"
                ).classes("font-semibold text-amber-900" + (" mt-1" if hidden else ""))
                for r in missing[:limit]:
                    ui.label(f"• {_row_text(r)}").classes("text-sm text-amber-900")
        if hidden and on_show_all is not None:
            ui.button("Show all candidates", on_click=on_show_all).props(
                "flat dense no-caps color=amber-10"
            ).mark("denovo-show-hidden")


def missing_validations(
    data_dir: Any, family_id: str, rows: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Validated present de novo calls of the family with no row in ``rows``."""
    curation = denovo.load_curation(data_dir, [family_id])
    keys = {
        (family_id, r.get("Variant") or "", str(r.get("sample") or "")) for r in rows
    }
    out = []
    for key, (status, inheritance) in sorted(curation.validations.items()):
        if key in keys or inheritance != "de novo":
            continue
        if status not in ("present", "in phase MNV"):
            continue
        out.append(
            {
                "Variant": key[1],
                "sample": key[2],
                "validation": status,
                "diagnostic": curation.diagnostics.get(key, ""),
            }
        )
    return out


def render_denovo_qc_badge(
    qc_rows: List[Dict[str, Any]], on_click: Optional[Callable[[], Any]] = None
) -> None:
    """Family header: a badge when a child's de novo count is out of range.

    Nothing when every child is ``ok`` or there is no verdict (family tables
    only): the warning is for what needs a look.
    """
    alerts = denovo.qc_alerts(qc_rows)
    if not alerts:
        return
    worst = "high" if any(v == "high" for _, v in alerts) else alerts[0][1]
    info = denovo.qc_info(worst)
    lines = []
    for child, verdict in alerts:
        q = next((r for r in qc_rows if str(r.get("child")) == child), {})
        rate = q.get("rate_per_gb")
        low, high = q.get("qc_low"), q.get("qc_high")
        text = f"{child}: {denovo.qc_info(verdict)['label']}"
        if isinstance(rate, (int, float)):
            text += f", {rate:.1f} strict per Gb"
        if isinstance(low, (int, float)) and isinstance(high, (int, float)):
            text += f" (range {low:.1f}–{high:.1f})"
        if q.get("sampleqc_flags"):
            text += f"; Sample QC: {q.get('sampleqc_flags')}"
        lines.append(text)
    badge = (
        ui.badge(f"De novo QC: {', '.join(c for c, _ in alerts)}", color=info["color"])
        .classes("cursor-pointer")
        .mark("denovo-qc-badge")
    )
    badge.tooltip(
        "Strict de novo count out of range (open SNVs/Indels):\n" + "\n".join(lines)
    ).props("max-width=40rem").classes("whitespace-pre-line")
    if on_click is not None:
        badge.on("click", lambda: on_click())
