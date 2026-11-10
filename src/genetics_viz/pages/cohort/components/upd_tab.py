"""UPD tab for the family page: uniparental disomy candidates, reviewed in one place.

One row per member and chromosome that any step points at (the upd step's
events, the ROH ``upd_candidate_chrom``, a Sample QC ``UPD_CANDIDATE`` flag)
or that a curator validated (:func:`genetics_viz.utils.upd.upd_candidates`).
The review dialog (:func:`genetics_viz.components.upd_dialog.show_upd_dialog`)
shows the evidence and the chromosome view, and records the validation and the
diagnostic. The tab has three blocks:

1. the candidates, with the upd step's newest label by default (a selector
   when several were run) and its High / Medium events (Low on demand);
2. the genome, each member's windows coloured by parent of origin (from
   ``upd_v4``, whose tables hold them);
3. the members: analysis mode and status, worded for the role.

Imprinting is the step's own ``imprinting`` column, as written.
"""

import logging
from typing import Any, Callable, Dict, List, Optional

from nicegui import ui

from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.offthread_panel import OffThreadPanel
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.upd_dialog import VALIDATION_COLORS, show_upd_dialog
from genetics_viz.utils import upd_calls
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.diagnostic_files import diagnostic_path
from genetics_viz.utils.export import export_filename, rows_to_tsv
from genetics_viz.utils.str_eh import diagnostic_badge
from genetics_viz.utils.upd import (
    UpdCandidate,
    family_upd_hints,
    read_upd_validations,
    upd_candidates,
)
from genetics_viz.utils.upd_calls import (
    CONFIDENCES,
    DEFAULT_CONFIDENCES,
    confidence_info,
    event_type,
    load_family_upd,
    load_windows,
    member_status_label,
    mode_info,
    origin_info,
    probe_upd_data,
    span_text,
    table_export,
    windows_of,
)
from genetics_viz.utils.upd_plots import genome_svg

logger = logging.getLogger(__name__)

__all__ = ["probe_upd_data", "render_upd_tab"]

_SOURCE_COLORS = {"ROH": "#6a1b9a", "Sample QC": "#f9a825", "Curator": "#757575"}


def _float(value: Any) -> Optional[float]:
    try:
        return None if value is None else round(float(value), 4)
    except (TypeError, ValueError):
        return None


def _validation_badge(cand: UpdCandidate) -> List[Dict[str, str]]:
    status = cand.status
    if not status:
        return [{"label": "not reviewed", "color": "#bdbdbd", "textColor": "#374151"}]
    origin = (cand.validation or {}).get("Inheritance") or ""
    label = status + (f" · {origin}" if origin and origin != "unknown" else "")
    tip = "Validated by {} on {}".format(
        (cand.validation or {}).get("User", ""),
        ((cand.validation or {}).get("Timestamp") or "")[:10],
    )
    return [
        {
            "label": label,
            "color": VALIDATION_COLORS.get(status, "#757575"),
            "tooltip": tip,
        }
    ]


def candidate_rows(
    candidates: List[UpdCandidate], diagnostics: Dict[Any, Any]
) -> List[Dict[str, Any]]:
    """The candidates table's rows, in the candidates' order."""
    rows = []
    for i, c in enumerate(candidates):
        e = c.event or {}
        info = confidence_info(e.get("confidence")) if e else None
        sources = []
        if e:
            sources.append(
                {
                    "label": f"UPD · {info['label']}",
                    "color": info["color"],
                    "tooltip": f"{len(c.events)} event(s) of the upd step",
                }
            )
        if c.roh:
            sources.append(
                {
                    "label": "ROH",
                    "color": _SOURCE_COLORS["ROH"],
                    "tooltip": "The ROH step's upd_candidate_chrom",
                }
            )
        if c.flags:
            sources.append(
                {
                    "label": "Sample QC",
                    "color": _SOURCE_COLORS["Sample QC"],
                    "textColor": "#111827",
                    "tooltip": "A UPD_CANDIDATE flag of Sample QC",
                }
            )
        if not sources:
            sources.append({"label": "Curator", "color": _SOURCE_COLORS["Curator"]})
        notes = []
        if e.get("recurrent") == "yes":
            notes.append({"label": "recurrent", "color": "#616161"})
        if e.get("shared_with"):
            notes.append(
                {
                    "label": f"shared with {e['shared_with']}",
                    "color": "#616161",
                    "tooltip": "A sibling has an overlapping event: an artefact",
                }
            )
        origin = origin_info(
            (c.validation or {}).get("Inheritance") or e.get("origin") or "unknown"
        )
        diag: Dict[str, Any] = {}
        add_diagnostic_status_to_row(diag, diagnostics, c.key, c.sample)
        diag_text = diag.get("Diagnostic") or ""
        span = c.span
        rows.append(
            {
                "key": f"{c.sample}:{c.chrom}",
                "_index": i,
                "sample": c.sample,
                "chrom": f"chr{c.chrom}",
                "_order": i,
                "sources": " ".join(s["label"] for s in sources),
                "sources_badges": sources,
                "span": span_text(*span) if span else "",
                "length_mb": _float(e.get("length_mb")),
                "type": event_type(e) if e else "",
                "origin": origin["label"],
                "origin_color": origin["color"],
                "imprinting": e.get("imprinting") or "",
                "lod": _float(e.get("lod")),
                "support": _float(e.get("support_fraction")),
                "depth": _float(e.get("depth_ratio")),
                "notes": " ".join(n["label"] for n in notes),
                "notes_badges": notes,
                "validation": c.status or "not reviewed",
                "validation_badges": _validation_badge(c),
                "diagnostic": diag_text,
                "diagnostic_badges": (
                    [
                        {
                            "label": diag_text,
                            "color": diagnostic_badge(diag_text)["color"],
                        }
                    ]
                    if diag_text
                    else []
                ),
            }
        )
    return rows


_COLUMNS: List[Dict[str, Any]] = [
    {
        "id": "actions",
        # Not "": the table falls back to the id for an empty header.
        "header": " ",
        "cellType": "action",
        "actionName": "review_upd",
        "actionIcon": "rate_review",
        "actionColor": "#1976d2",
        "actionTooltip": "Evidence, chromosome view, validation and diagnostic",
        "sortable": False,
    },
    {"id": "sample", "header": "Member", "sortable": True},
    {"id": "chrom", "header": "Chromosome", "sortable": True, "sortField": "_order"},
    {"id": "sources", "header": "Found by", "cellType": "badge_list"},
    {"id": "span", "header": "Span", "sortable": True},
    {"id": "type", "header": "Type"},
    {"id": "origin", "header": "Origin", "cellType": "badge", "colorField": "origin_color"},
    {"id": "imprinting", "header": "Imprinting", "maxWidth": 260},
    {"id": "lod", "header": "LOD", "cellType": "number", "sortable": True},
    {"id": "support", "header": "Support", "cellType": "number", "sortable": True},
    {"id": "depth", "header": "Depth", "cellType": "number", "sortable": True},
    {"id": "notes", "header": "Notes", "cellType": "badge_list"},
    {"id": "validation", "header": "Validation", "cellType": "badge_list"},
    {"id": "diagnostic", "header": "Diagnostic", "cellType": "badge_list"},
]  # fmt: skip


def _section(title: str, help_text: str) -> None:
    with ui.row().classes("items-center gap-2 mt-4"):
        ui.label(title).classes("text-base font-semibold text-blue-700")
        ui.icon("info_outline", color="grey").classes("text-sm cursor-help").tooltip(
            help_text
        )


def _render_members(
    calls: upd_calls.UpdData, members: List[str], show_detail: bool
) -> None:
    rows = [m for m in calls.members if m.get("sample") in members]
    if not rows:
        return
    with (
        ui.expansion("Members analysed", icon="groups")
        .classes("w-full mt-2")
        .mark("upd-members")
    ):
        columns = [
            ("sample", "Member"),
            ("mode", "Mode"),
            *([("parents", "Parents used")] if show_detail else []),
            ("chrx", "chrX analysed"),
            ("n_sites", "Sites"),
            ("error_rate", "Error rate"),
            ("status", "Status"),
            ("n_events", "Events"),
        ]
        ui.table(
            columns=[
                {"name": k, "label": label, "field": k, "align": "left"}
                for k, label in columns
            ],
            rows=[
                {
                    "sample": m.get("sample"),
                    "mode": mode_info(m.get("mode"))["label"],
                    "parents": " + ".join(
                        p for p in (m.get("father_used"), m.get("mother_used")) if p
                    )
                    or "—",
                    "chrx": m.get("chrx") or "",
                    "n_sites": f"{m['n_sites']:,}" if m.get("n_sites") else "—",
                    "error_rate": m.get("error_rate") if m.get("error_rate") else "—",
                    "status": member_status_label(m.get("status"), show_detail),
                    "n_events": m.get("n_events") or 0,
                }
                for m in rows
            ],
            row_key="sample",
        ).props("dense flat").classes("w-full")
        ui.label(
            " ".join(
                f"{mode_info(k)['label']}: {mode_info(k)['help']}"
                for k in ("trio", "duo", "solo")
            )
        ).classes("text-xs text-gray-500")


def _render_genome(
    calls: upd_calls.UpdData, members: List[str], confidences: List[str]
) -> None:
    _section(
        "Genome",
        "Each member's 500 kb windows: pink where the child's origin-informative"
        " alleles come from the mother only, blue from the father only, grey when"
        " balanced; dark purple where a member without parents has no het call."
        " Events are outlined by confidence.",
    )
    panel = OffThreadPanel("upd-genome")

    def load() -> Dict[str, Dict[str, List[Dict[str, Any]]]]:
        windows = load_windows(calls)
        out: Dict[str, Dict[str, List[Dict[str, Any]]]] = {}
        for m in members:
            by_chrom: Dict[str, List[Dict[str, Any]]] = {}
            for w in windows_of(windows, m):
                by_chrom.setdefault(w["chrom"], []).append(w)
            if by_chrom:
                out[m] = by_chrom
        return out

    def render(by_member: Dict[str, Dict[str, List[Dict[str, Any]]]]) -> None:
        if not by_member:
            ui.label("No windows for the selected members.").classes(
                "text-gray-500 italic text-sm"
            )
            return
        lanes = [m for m in members if m in by_member]
        modes = {
            m: mode_info((calls.member(m) or {}).get("mode"))["label"] for m in lanes
        }
        events = {
            m: [e for e in calls.events_of(m) if e.get("confidence") in confidences]
            for m in lanes
        }
        ui.label("Lanes: " + ", ".join(f"{m} ({modes[m]})" for m in lanes)).classes(
            "text-xs text-gray-500"
        )
        ui.html(genome_svg(lanes, by_member, events, modes), sanitize=False).classes(
            "w-full"
        ).mark("upd-genome-svg")

    panel.load(load, render, "Loading the genome windows…")


def render_upd_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
    on_acks_changed: Optional[Callable[[], None]] = None,
    link: Any = None,
    upd_listeners: Optional[List[Callable[[], None]]] = None,
) -> Callable[[str, str], None]:
    """Render the UPD tab panel content; returns ``open_review(sample, chrom)``.

    Args:
        on_acks_changed: Called after a UPD validation or diagnostic, which
            may acknowledge a Sample QC flag, so the page header follows.
        link: The page's :class:`PageLink` (``updl`` label, ``updc``
            confidences, ``upd`` the review open).
        upd_listeners: Called after a UPD review anywhere on the page; the
            tab adds its own redraw.
    """
    data_dir = store.data_dir
    show_detail = can_write()
    linked_label = link.get("updl") if link is not None else None
    calls = load_family_upd(data_dir, family_id, linked_label)
    linked_conf = link.get_list("updc") if link is not None else None
    state: Dict[str, Any] = {
        "calls": calls,
        "confidences": [
            c for c in (linked_conf or DEFAULT_CONFIDENCES) if c in CONFIDENCES
        ],
        "hints": family_upd_hints(data_dir, family_id),
    }

    def push_link() -> None:
        if link is None:
            return
        calls_now = state["calls"]
        link.set(
            "updl",
            calls_now.label if calls_now is not None else None,
            default=calls_now.labels[0] if calls_now is not None else None,
        )
        link.set("updc", state["confidences"], default=DEFAULT_CONFIDENCES)
        link.push()

    with ui.row().classes("items-center gap-2 mb-2 w-full"):
        ui.label("Uniparental disomy").classes("text-lg font-semibold text-blue-700")
        provenance = ui.label(
            f"upd step · {calls.label}" if calls is not None else "no upd step results"
        ).classes("text-xs text-gray-500")
        ui.space()

        def members() -> List[str]:
            return list(selected_members["value"])

        def build(kind: str) -> Callable[[], str]:
            def builder() -> str:
                current = state["calls"]
                if current is None:
                    return ""
                wanted = set(members())
                if kind == "windows":
                    df = load_windows(current)
                    rows = [r for m in members() for r in windows_of(df, m)]
                else:
                    rows = [
                        r for r in getattr(current, kind) if r.get("sample") in wanted
                    ]
                out, columns = table_export(rows, kind, show_detail)
                return rows_to_tsv(out, columns)

            return builder

        def filename(kind: Optional[str], ext: str = "tsv") -> Callable[[], str]:
            return lambda: export_filename(
                family_id,
                "upd",
                kind,
                tag=state["calls"].label if state["calls"] is not None else None,
                ext=ext,
            )

        if calls is not None:
            downloads = [
                Download("events", "UPD events", filename("events"), build("events")),
                Download(
                    "segments", "UPD segments", filename("segments"), build("segments")
                ),
                Download("members", "Members", filename("members"), build("members")),
            ]
            if calls.has_windows:
                downloads.append(
                    Download(
                        "windows",
                        "Windows (genotype counts)",
                        filename("windows"),
                        build("windows"),
                    )
                )
            download_button(
                downloads, mark="download-upd", zip_filename=filename(None, "zip")
            )

    with ui.row().classes("items-center gap-3 mb-2 w-full"):
        if calls is not None and len(calls.labels) > 1:

            def on_label(event: Any) -> None:
                state["calls"] = load_family_upd(data_dir, family_id, event.value)
                provenance.text = f"upd step · {state['calls'].label}"
                push_link()
                render_content.refresh()

            ui.select(
                options={
                    label: label + (" (newest)" if i == 0 else "")
                    for i, label in enumerate(calls.labels)
                },
                value=calls.label,
                label="upd version",
                on_change=on_label,
            ).props("outlined dense").classes("w-56").mark("upd-label")

        if calls is not None:
            ui.label("Events:").classes("text-sm text-gray-600")
            for key in CONFIDENCES:

                def make_toggle(k: str) -> Callable[[Any], None]:
                    def handler(e: Any) -> None:
                        chosen = set(state["confidences"])
                        if e.value:
                            chosen.add(k)
                        else:
                            chosen.discard(k)
                        state["confidences"] = [c for c in CONFIDENCES if c in chosen]
                        push_link()
                        render_content.refresh()

                    return handler

                ui.checkbox(
                    confidence_info(key)["label"],
                    value=key in state["confidences"],
                    on_change=make_toggle(key),
                ).props("dense").mark(f"upd-confidence-{key}")

    if calls is None:
        ui.label(
            "The upd step has not run for this family. The candidates below come"
            " from the ROH and Sample QC steps."
        ).classes("text-sm text-gray-500 mb-2")

    # Outside the refreshable section: a save redraws it, which would otherwise
    # delete the review dialog.
    review_anchor = ui.element("div")

    def refresh() -> None:
        state["hints"] = family_upd_hints(data_dir, family_id)
        render_content.refresh()

    def saved() -> None:
        # The page's on_acks_changed runs every UPD listener, this tab's too.
        if on_acks_changed is None or upd_listeners is None:
            refresh()
        if on_acks_changed is not None:
            on_acks_changed()

    def closed() -> None:
        if link is not None:
            link.set("upd", None)
            link.push()

    def open_review(sample: str, chrom: str) -> None:
        if link is not None:
            link.set("upd", f"{sample}:{upd_calls.normalise_chrom(chrom)}")
            link.push()
        current = state["calls"]
        show_upd_dialog(
            review_anchor,
            data_dir,
            family_id,
            sample,
            chrom,
            on_saved=saved,
            on_closed=closed,
            upd_label=current.label if current is not None else None,
        )

    @ui.refreshable
    def render_content() -> None:
        current: Optional[upd_calls.UpdData] = state["calls"]
        member_ids = members()
        if not member_ids:
            ui.label("No members selected").classes("text-gray-500 italic")
            return
        hints = state["hints"]
        candidates = upd_candidates(
            current.events if current is not None else [],
            hints.froh_rows,
            hints.flags,
            read_upd_validations(data_dir, family_id),
            state["confidences"],
            member_ids,
        )
        _section(
            "Candidates",
            "One row per member and chromosome that the upd step (at the ticked"
            " confidences), the ROH step or Sample QC points at, or that a"
            " curator validated. Review a row for its evidence, the chromosome"
            " view, the validation and the diagnostic.",
        )
        if not candidates:
            ui.label(
                "No UPD candidate among the selected members"
                + (" at these confidences." if current is not None else ".")
            ).classes("text-gray-500 italic text-sm").mark("upd-no-candidates")
        else:
            diagnostics = load_diagnostic_map(
                diagnostic_path(data_dir, "UPD"), family_id
            )
            rows = candidate_rows(candidates, diagnostics)

            def on_row_action(e: Dict[str, Any]) -> None:
                index = (e.get("row") or {}).get("_index")
                if isinstance(index, int) and 0 <= index < len(candidates):
                    cand = candidates[index]
                    open_review(cand.sample, cand.chrom)

            with ui.element("div").classes("w-full").mark("upd-candidates"):
                DataTable(
                    columns=_COLUMNS,
                    rows=rows,
                    row_key="key",
                    pagination={"rowsPerPage": 20},
                    on_row_action=on_row_action,
                )
            ui.label(
                "A validation (present or absent) acknowledges the matching Sample"
                " QC flag and settles the ROH candidate. Imprinting is as the upd"
                " step words it, per chromosome or region."
            ).classes("text-xs text-gray-500")

        if current is not None and current.has_windows:
            _render_genome(current, member_ids, state["confidences"])
        elif current is not None:
            ui.label(
                f"{current.label} wrote no per-window genotype counts (from upd_v4):"
                " the genome overview and the chromosome view's bars need them."
            ).classes("text-xs text-gray-500 mt-2")
        if current is not None:
            _render_members(current, member_ids, show_detail)

    render_content()
    data_table_refreshers.append(render_content.refresh)
    if upd_listeners is not None:
        upd_listeners.append(refresh)
    return open_review
