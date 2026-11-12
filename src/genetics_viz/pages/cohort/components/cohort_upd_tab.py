"""UPD sub-tab of the cohort page's right-hand panel, and the individuals table's UPD column.

The cohort's upd tables (``cohorts/<NAME>/upd/<NAME>.<label>.*.tsv``, the
family tables concatenated, with ``recurrent`` on the events), restricted to
the individuals the table's filters leave:

- what was analysed, by mode, and how many events at each confidence;
- events per chromosome, stacked by confidence (recurrent regions stand out);
- the events, each linking to its family's UPD tab with the review open, with
  their standing validation.

The newest label is shown by default; the selector lists the others.
"""

import asyncio
import logging
from typing import Any, Callable, Dict, List, Optional, Set, Tuple
from urllib.parse import quote, urlencode

from nicegui import ui

from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.upd_dialog import VALIDATION_COLORS
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.export import export_filename, rows_to_tsv
from genetics_viz.utils.upd import read_upd_validations, standing_validations
from genetics_viz.utils.upd_calls import (
    CONFIDENCES,
    DEFAULT_CONFIDENCES,
    UpdData,
    confidence_info,
    event_type,
    load_cohort_upd,
    mode_info,
    origin_info,
    probe_cohort_upd_data,
    span_text,
    table_export,
)
from genetics_viz.utils.upd_plots import CHROMOSOMES, cohort_chromosome_option

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_upd_data", "render_cohort_upd_tab", "upd_column_fields"]

Key = Tuple[str, str]
#: ``(FID, sample, chrom) -> standing validation``.
Standing = Dict[Tuple[str, str, str], Dict[str, str]]

#: The individuals table's UPD filter options besides the confidences.
UPD_ANY = "Any event"
UPD_PRESENT = "Validated present"
UPD_ABSENT = "Validated absent"
UPD_TO_REVIEW = "Not reviewed"
UPD_NONE = "No event"
UPD_NOT_ANALYSED = "Not analysed"


def upd_filter_options() -> List[str]:
    return [
        UPD_ANY,
        *(confidence_info(c)["label"] for c in CONFIDENCES),
        UPD_PRESENT,
        UPD_ABSENT,
        UPD_TO_REVIEW,
        UPD_NONE,
        UPD_NOT_ANALYSED,
    ]


def cohort_standing(data_dir: Any) -> Standing:
    """Every family's standing UPD validations."""
    rows = read_upd_validations(data_dir)
    by_family: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        by_family.setdefault(row.get("FID") or "", []).append(row)
    out: Standing = {}
    for fid, family_rows in by_family.items():
        for (sample, chrom), row in standing_validations(family_rows).items():
            out[(fid, sample, chrom)] = row
    return out


def upd_link(cohort_name: str, fid: str, sample: str, chrom: str) -> str:
    """The family page's UPD tab with the review of ``sample``'s ``chrom`` open."""
    query = urlencode({"tab": "upd", "upd": f"{sample}:{chrom}"})
    return f"/cohort/{quote(cohort_name)}/family/{quote(fid)}?{query}"


def upd_column_fields(
    calls: Optional[UpdData], standing: Standing
) -> Callable[[str, str], Dict[str, Any]]:
    """The individuals table's ``UPD`` cell of each ``(FID, sample)``.

    One badge per chromosome with an event of the default confidences: the
    validation's colour once a curator settled it, the confidence's before.
    ``_upd_open`` lists the confidence of each such chromosome no curator
    validated present or absent: the UPD tab's badge counts them.
    """
    if calls is None:
        return lambda fid, sid: {}
    analysed = {
        (m.get("family_id") or "", m.get("sample") or "") for m in calls.members
    }
    by_key: Dict[Key, Dict[str, Dict[str, Any]]] = {}
    for e in calls.events:
        key = (e.get("family_id") or "", e.get("sample") or "")
        best = by_key.setdefault(key, {}).get(e.get("chrom") or "")
        if (
            best is None
            or confidence_info(e.get("confidence"))["rank"]
            > confidence_info(best.get("confidence"))["rank"]
        ):
            by_key[key][e.get("chrom") or ""] = e
    order = {c: i for i, c in enumerate(CHROMOSOMES)}

    def fields(fid: str, sid: str) -> Dict[str, Any]:
        if (fid, sid) not in analysed and (fid, sid) not in by_key:
            return {
                "UPD": "",
                "_upd_cats": [UPD_NOT_ANALYSED],
                "_upd_rank": -1,
                "_upd_open": [],
            }
        badges, cats, rank, open_ = [], set(), 0, []
        for chrom, e in sorted(
            by_key.get((fid, sid), {}).items(), key=lambda t: order.get(t[0], 99)
        ):
            info = confidence_info(e.get("confidence"))
            cats.add(info["label"])
            row = standing.get((fid, sid, chrom))
            status = (row or {}).get("Validation", "")
            if status == "present":
                cats.add(UPD_PRESENT)
            elif status == "absent":
                cats.add(UPD_ABSENT)
            else:
                cats.add(UPD_TO_REVIEW)
            if e.get("confidence") not in DEFAULT_CONFIDENCES and not status:
                continue
            if status not in ("present", "absent"):
                open_.append(e.get("confidence") or "")
            rank = max(rank, info["rank"])
            label = f"chr{chrom} · {status or info['label']}"
            badges.append(
                {
                    "label": label,
                    "color": VALIDATION_COLORS.get(status, info["color"])
                    if status
                    else info["color"],
                    "tooltip": (
                        f"{span_text(e.get('start'), e.get('end'))},"
                        f" {origin_info(e.get('origin'))['label']}"
                        f" {event_type(e)} ({info['label']})"
                        + (f"; validated {status}" if status else "; not reviewed")
                    ),
                }
            )
        if by_key.get((fid, sid)):
            cats.add(UPD_ANY)
        else:
            cats.add(UPD_NONE)
        return {
            "UPD": "",
            "UPD_badges": badges,
            "_upd_cats": sorted(cats),
            "_upd_rank": rank,
            "_upd_open": open_,
        }

    return fields


def _keys(filtered_state: Dict[str, Any]) -> Set[Key]:
    return {
        (ind.get("FID") or "", ind.get("Sample ID") or "")
        for ind in filtered_state["individuals"]
        if ind.get("Sample ID")
    }


def _event_rows(
    events: List[Dict[str, Any]], standing: Standing, cohort_name: str
) -> List[Dict[str, Any]]:
    order = {c: i for i, c in enumerate(CHROMOSOMES)}
    rows = []
    for i, e in enumerate(events):
        fid, sample, chrom = (
            e.get("family_id") or "",
            e.get("sample") or "",
            e.get("chrom") or "",
        )
        info = confidence_info(e.get("confidence"))
        origin = origin_info(e.get("origin"))
        row = standing.get((fid, sample, chrom))
        status = (row or {}).get("Validation", "")
        rows.append(
            {
                "key": f"{fid}|{sample}|{e.get('event_id') or i}",
                "FID": fid,
                "_href": upd_link(cohort_name, fid, sample, chrom),
                "sample": sample,
                "chrom": f"chr{chrom}",
                "_chrom_order": order.get(chrom, 99) * 1e10 + (e.get("start") or 0),
                "span": span_text(e.get("start"), e.get("end")),
                "length_mb": e.get("length_mb"),
                "origin": origin["label"],
                "origin_color": origin["color"],
                "type": event_type(e),
                "mode": mode_info(e.get("mode"))["label"],
                "confidence": info["label"],
                "confidence_color": info["color"],
                "_rank": info["rank"],
                "imprinting": e.get("imprinting") or "",
                "notes": "",
                "notes_badges": (
                    [
                        {
                            "label": "recurrent",
                            "color": "#616161",
                            "tooltip": f"Overlaps events of {e.get('recurrent_families')}"
                            " other families: most likely an artefact region",
                        }
                    ]
                    if e.get("recurrent") == "yes"
                    else []
                ),
                "validation": status or "not reviewed",
                "validation_badges": [
                    {
                        "label": status
                        + (
                            f" · {row.get('Inheritance')}"
                            if row and row.get("Inheritance") not in ("", "unknown")
                            else ""
                        ),
                        "color": VALIDATION_COLORS.get(status, "#757575"),
                    }
                ]
                if status
                else [
                    {
                        "label": "not reviewed",
                        "color": "#bdbdbd",
                        "textColor": "#374151",
                    }
                ],
            }
        )
    rows.sort(key=lambda r: (-r["_rank"], r["_chrom_order"], r["FID"], r["sample"]))
    return rows


def _render_summary(
    calls: UpdData, keys: Set[Key], events: List[Dict[str, Any]], standing: Standing
) -> None:
    members = [
        m
        for m in calls.members
        if (m.get("family_id") or "", m.get("sample") or "") in keys
    ]
    modes: Dict[str, int] = {}
    for m in members:
        label = mode_info(m.get("mode"))["label"]
        modes[label] = modes.get(label, 0) + 1
    by_conf: Dict[str, int] = {}
    for e in events:
        by_conf[e.get("confidence") or ""] = (
            by_conf.get(e.get("confidence") or "", 0) + 1
        )
    chromosomes = {
        (e.get("family_id") or "", e.get("sample") or "", e.get("chrom") or "")
        for e in events
    }
    reviewed = sum(
        1
        for k in chromosomes
        if (standing.get(k) or {}).get("Validation") in ("present", "absent")
    )
    with ui.column().classes("w-full gap-1 text-sm").mark("cohort-upd-summary"):
        ui.label(
            f"{len(members)} individuals analysed"
            + (
                " ("
                + ", ".join(f"{n} {mode}" for mode, n in sorted(modes.items()))
                + ")"
                if modes
                else ""
            )
        )
        with ui.row().classes("items-center gap-2"):
            ui.label(f"{len(events)} events:")
            for conf in CONFIDENCES:
                n = by_conf.get(conf, 0)
                if n:
                    info = confidence_info(conf)
                    ui.label(f"{n} {info['label'].lower()}").classes(
                        "text-xs text-white rounded px-1"
                    ).style(f"background-color: {info['color']}")
        ui.label(
            f"{reviewed} of {len(chromosomes)} candidate chromosomes reviewed"
            " (validated present or absent)"
        ).classes("text-xs text-gray-600")


_COLUMNS: List[Dict[str, Any]] = [
    {"id": "sample", "header": "Sample", "sortable": True},
    {
        "id": "FID",
        "header": "Family",
        "cellType": "link",
        "hrefField": "_href",
        "sortable": True,
    },
    {"id": "chrom", "header": "Chr", "sortable": True, "sortField": "_chrom_order",
     "sorting": "numerical"},
    {"id": "span", "header": "Span"},
    {"id": "confidence", "header": "Confidence", "cellType": "badge",
     "colorField": "confidence_color", "sortable": True, "sortField": "_rank",
     "sorting": "numerical"},
    {"id": "origin", "header": "Origin", "cellType": "badge", "colorField": "origin_color"},
    {"id": "type", "header": "Type", "maxWidth": 140},
    {"id": "mode", "header": "Mode", "sortable": True},
    # Long, and the panel is narrow: wrapped.
    {"id": "imprinting", "header": "Imprinting", "maxWidth": 200},
    {"id": "notes", "header": "Notes", "cellType": "badge_list"},
    {"id": "validation", "header": "Validation", "cellType": "badge_list", "sortable": True},
]  # fmt: skip


def render_cohort_upd_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort UPD panel; it follows the individuals table's filters."""
    show_detail = can_write()
    state: Dict[str, Any] = {
        "calls": None,
        "label": None,
        "confidences": list(DEFAULT_CONFIDENCES),
        "recurrent": True,
        "loading": True,
        "standing": {},
    }

    def shown_events() -> List[Dict[str, Any]]:
        calls: Optional[UpdData] = state["calls"]
        if calls is None:
            return []
        keys = _keys(filtered_state)
        return [
            e
            for e in calls.events
            if (e.get("family_id") or "", e.get("sample") or "") in keys
            and e.get("confidence") in state["confidences"]
            and (state["recurrent"] or e.get("recurrent") != "yes")
        ]

    def build(kind: str) -> Callable[[], str]:
        def builder() -> str:
            calls: Optional[UpdData] = state["calls"]
            if calls is None:
                return ""
            if kind == "events":
                rows = shown_events()
            else:
                keys = _keys(filtered_state)
                rows = [
                    r
                    for r in getattr(calls, kind)
                    if (r.get("family_id") or "", r.get("sample") or "") in keys
                ]
            out, columns = table_export(rows, kind, show_detail)
            return rows_to_tsv(out, columns)

        return builder

    def filename(kind: Optional[str], ext: str = "tsv") -> Callable[[], str]:
        return lambda: export_filename(
            cohort.name, "upd", kind, tag=state["label"], ext=ext
        )

    with ui.card().classes("w-full"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("Uniparental disomy").classes(
                "text-lg font-semibold text-blue-700"
            )
            provenance = ui.label("").classes("text-xs text-gray-500")
            ui.space()
            download_button(
                [
                    Download("events", "UPD events (shown)", filename("events"), build("events")),
                    Download("segments", "UPD segments", filename("segments"), build("segments")),
                    Download("members", "Members", filename("members"), build("members")),
                ],
                mark="download-upd",
                zip_filename=filename(None, "zip"),
            )  # fmt: skip

        @ui.refreshable
        def controls() -> None:
            calls: Optional[UpdData] = state["calls"]
            if calls is None:
                return
            with ui.row().classes("items-center gap-3 w-full flex-wrap"):
                if len(calls.labels) > 1:
                    ui.select(
                        options={
                            label: label + (" (newest)" if i == 0 else "")
                            for i, label in enumerate(calls.labels)
                        },
                        value=calls.label,
                        label="upd version",
                        on_change=lambda e: reload(e.value),
                    ).props("outlined dense").classes("w-48").mark("cohort-upd-label")
                for key in CONFIDENCES:

                    def make_toggle(k: str) -> Callable[[Any], None]:
                        def handler(e: Any) -> None:
                            chosen = set(state["confidences"])
                            if e.value:
                                chosen.add(k)
                            else:
                                chosen.discard(k)
                            state["confidences"] = [
                                c for c in CONFIDENCES if c in chosen
                            ]
                            render_body.refresh()

                        return handler

                    ui.checkbox(
                        confidence_info(key)["label"],
                        value=key in state["confidences"],
                        on_change=make_toggle(key),
                    ).props("dense").mark(f"cohort-upd-confidence-{key}")

                def on_recurrent(e: Any) -> None:
                    state["recurrent"] = bool(e.value)
                    render_body.refresh()

                ui.switch(
                    "Include recurrent",
                    value=state["recurrent"],
                    on_change=on_recurrent,
                ).props("dense").tooltip(
                    "Events overlapping those of three or more other families: most"
                    " likely artefact regions"
                ).mark("cohort-upd-recurrent")

        controls()

        @ui.refreshable
        def render_body() -> None:
            if state["loading"]:
                with ui.row().classes("items-center gap-2 p-4"):
                    ui.spinner(size="md")
                    ui.label("Loading UPD calls…").classes("text-gray-500")
                return
            calls: Optional[UpdData] = state["calls"]
            if calls is None:
                ui.label("No cohort-level UPD results").classes("text-gray-500 italic")
                return
            keys = _keys(filtered_state)
            events = shown_events()
            _render_summary(calls, keys, events, state["standing"])
            if not events:
                ui.label(
                    "No event among the filtered individuals at these confidences."
                ).classes("text-gray-500 italic text-sm")
                return
            ui.echart(cohort_chromosome_option(events)).classes("w-full").style(
                "height: 220px"
            ).mark("cohort-upd-chart")
            with ui.element("div").classes("w-full").mark("cohort-upd-events"):
                DataTable(
                    columns=_COLUMNS,
                    rows=_event_rows(events, state["standing"], cohort.name),
                    row_key="key",
                    pagination={"rowsPerPage": 10},
                )
            ui.label(
                "A family opens on its UPD tab with the review of that chromosome."
                " Imprinting is as the upd step words it."
            ).classes("text-xs text-gray-500")

        render_body()

    async def load(label: Optional[str] = None) -> None:
        try:
            calls = await asyncio.to_thread(load_cohort_upd, cohort, label)
            standing = await asyncio.to_thread(cohort_standing, store.data_dir)
        except Exception:
            logger.warning(
                "Failed to read the UPD calls of %s", cohort.name, exc_info=True
            )
            calls, standing = None, {}
        state.update(calls=calls, standing=standing, loading=False)
        state["label"] = calls.label if calls is not None else None
        provenance.text = f"upd step · {calls.label}" if calls is not None else ""
        controls.refresh()
        render_body.refresh()

    def reload(label: str) -> None:
        state["loading"] = True
        render_body.refresh()

        async def run() -> None:
            await load(label)

        ui.timer(0.05, run, once=True)

    refreshers.append(render_body.refresh)
    ui.timer(0.05, load, once=True)
