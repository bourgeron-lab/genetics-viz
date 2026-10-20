"""Cohort-wide Sample QC sub-tab for the cohort page's right-hand panel.

Follows the filters of the individuals table:

1. how many of the shown individuals carry an issue, per category;
2. the flagged individuals, linked to their family;
3. relatives across families: a count, the hubs, and a dialog with the
   network of every cross-family related pair touching a shown individual;
4. the X/Y copy-number plane, the KING / IBS0 plane, contamination and
   Mendelian errors, the flagged individuals marked.

The cohort tables are tens of KB; they are still read off-thread, as the
other cohort tabs' are. Wording follows the user's role (see ``sample_qc``).

Flags a curator acknowledged (``sample_qc_acks``) are left out of the counts,
the flagged individuals and the plot marks, unless "Include acknowledged" is
on; they are then listed greyed.
"""

import asyncio
import logging
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from nicegui import ui

from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.data import get_static_prefix
from genetics_viz.utils.export import (
    export_filename,
    people_from_individuals,
    rows_to_tsv,
)
from genetics_viz.utils.roh import to_float
from genetics_viz.utils.sample_qc import (
    HUB_MINIMUM,
    CATEGORIES,
    QUALITY_CATEGORIES,
    SEVERITIES,
    SEX_CATEGORIES,
    THRESHOLDS,
    FlagRow,
    PairRow,
    QcData,
    cohort_report_url,
    cross_family_graph,
    find_qc_labels,
    flag_badges,
    flags_by_key,
    get_cohort_qc_dir,
    load_cohort_qc,
    mendelian_fence,
    pair_rows,
    probe_cohort_qc_data,
    relationship_info,
    sample_point,
    severity_info,
)
from genetics_viz.utils.sample_qc_acks import (
    Ack,
    AckKey,
    acknowledged_badges,
    load_acks,
    split_acknowledged,
)
from genetics_viz.utils.sample_qc_export import (
    flags_export,
    pairs_export,
    samples_export,
)
from genetics_viz.utils.sample_qc_plots import (
    cross_family_graph_option,
    metric_strip_option,
    relatedness_scatter_option,
    sex_scatter_option,
)

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_qc_data", "render_cohort_sample_qc_tab"]

_NUM = {"cellType": "number", "sorting": "numerical", "sortable": True}


Key = Tuple[str, str]


def _table_box() -> ui.element:
    """A block box the table can scroll inside (see cohort_roh_tab)."""
    return ui.element("div").classes("w-full min-w-0")


def _keys(filtered_state: Dict[str, Any]) -> Set[Key]:
    return {
        (ind.get("FID") or "", ind["Sample ID"])
        for ind in filtered_state.get("individuals", [])
        if ind.get("Sample ID")
    }


def _subtitle(text: str) -> None:
    ui.label(text).classes("text-sm font-semibold text-blue-700 mt-3")


# ---------------------------------------------------------------------------
# Blocks
# ---------------------------------------------------------------------------


def _render_counts(by_key: Dict[Key, List[FlagRow]], keys: Set[Key]) -> None:
    """Individuals per category, by their worst severity in that category."""
    counts: Dict[str, Dict[str, int]] = {c: {} for c in CATEGORIES}
    for key in keys:
        worst: Dict[str, FlagRow] = {}
        for f in by_key.get(key, []):
            cat = f.info.category
            if cat not in worst or f.rank > worst[cat].rank:
                worst[cat] = f
        for cat, f in worst.items():
            sev = severity_info(f.severity)["key"]
            counts.setdefault(cat, {})[sev] = counts.setdefault(cat, {}).get(sev, 0) + 1
    severities = list(SEVERITIES)
    # The category names need the room; the counts are short.
    with (
        ui.element("div")
        .classes("w-full gap-x-3 gap-y-0 text-sm")
        .style(
            "display: grid; grid-template-columns:"
            f" minmax(0, 2.5fr) repeat({len(severities)}, minmax(0, 1fr))"
        )
    ):
        ui.label("Individuals with").classes("font-semibold text-gray-600")
        for sev in severities:
            ui.label(severity_info(sev)["label"]).classes("font-semibold").style(
                f"color: {severity_info(sev)['color']}"
            )
        for cat, cfg in CATEGORIES.items():
            ui.label(cfg.get("label", cat)).classes("truncate")
            for sev in severities:
                n = counts.get(cat, {}).get(sev, 0)
                ui.label(str(n) if n else "·").classes(
                    "font-semibold" if n else "text-gray-400"
                )


def _render_flagged(
    by_key: Dict[Key, List[FlagRow]],
    keys: Set[Key],
    cross_counts: Dict[Key, int],
    cohort_name: str,
    show_detail: bool,
    acks: Dict[AckKey, Ack],
) -> None:
    rows = []
    for key in keys:
        flags = by_key.get(key)
        if not flags:
            continue
        active, acked = split_acknowledged(flags, acks)
        # An individual with acknowledged flags only sorts below the others.
        worst = max(active, key=lambda f: f.rank).rank if active else 0
        rows.append(
            {
                "key": f"{key[0]}|{key[1]}",
                "FID": key[0],
                "sample": key[1],
                "issues": "",
                "issues_badges": flag_badges(active, show_detail)
                + acknowledged_badges(acked, acks, show_detail),
                "rank": worst * 1000 + len(active),
                "relatives": cross_counts.get(key, 0) or None,
            }
        )
    # Most severe first, matching the sort indicator.
    rows.sort(key=lambda r: (-r["rank"], r["FID"], r["sample"]))
    ui.label(f"{len(rows)} of {len(keys)} individuals flagged").classes(
        "text-xs text-gray-500"
    )
    if not rows:
        return
    with _table_box():
        DataTable(
            columns=[
                {"id": "sample", "header": "Sample", "sortable": True},
                {
                    "id": "FID",
                    "header": "Family",
                    "cellType": "link",
                    "href": f"/cohort/{cohort_name}/family/{{FID}}",
                    "sortable": True,
                },
                {
                    "id": "issues",
                    "header": "Issues",
                    "cellType": "badge_list",
                    "sortable": True,
                    "sortField": "rank",
                    "sorting": "numerical",
                },
                {"id": "relatives", "header": "Rel. in other fam.", **_NUM},
            ],
            rows=rows,
            row_key="key",
            pagination={"rowsPerPage": 10},
            initial_sorting=[{"id": "issues", "desc": True}],
        )


def _cross_pair_table(
    pairs: List[PairRow], cohort_name: str, show_detail: bool
) -> None:
    rows = [
        {
            "key": p.key,
            "sample_a": p.sample_a,
            "family_a": p.family_a,
            "sample_b": p.sample_b,
            "family_b": p.family_b,
            "observed": relationship_info(p.observed)["label"],
            "observed_color": relationship_info(p.observed)["color"],
            "kinship": p.kinship,
            "ratio": p.ratio,
            "n_sites": p.n_sites,
        }
        for p in sorted(pairs, key=lambda p: -(p.kinship or 0))
    ]
    link = f"/cohort/{cohort_name}/family/"
    columns: List[Dict[str, Any]] = [
        {"id": "sample_a", "header": "Sample", "sortable": True},
        {
            "id": "family_a",
            "header": "Family",
            "cellType": "link",
            "href": link + "{family_a}",
            "sortable": True,
        },
        {"id": "sample_b", "header": "Relative", "sortable": True},
        {
            "id": "family_b",
            "header": "Their family",
            "cellType": "link",
            "href": link + "{family_b}",
            "sortable": True,
        },
    ]
    if show_detail:
        columns.append(
            {
                "id": "observed",
                "header": "Observed",
                "cellType": "badge",
                "colorField": "observed_color",
                "sortable": True,
            }
        )
    columns += [
        {"id": "kinship", "header": "KING φ", **_NUM},
        {"id": "ratio", "header": "IBS0/IBS2", **_NUM},
        {"id": "n_sites", "header": "Sites", **_NUM},
    ]
    with _table_box():
        DataTable(
            columns=columns,
            rows=rows,
            row_key="key",
            pagination={"rowsPerPage": 15},
            initial_sorting=[{"id": "kinship", "desc": True}],
        )


def _open_cross_family_dialog(
    pairs: List[PairRow],
    keys: Set[Key],
    cohort_name: str,
    show_detail: bool,
) -> None:
    graph = cross_family_graph(pairs, restrict_to=keys)
    hubs = graph.hubs(HUB_MINIMUM)
    with (
        ui.dialog().props("maximized=false") as dialog,
        ui.card().classes("w-[95vw] max-w-none"),
    ):
        with ui.row().classes("w-full items-center"):
            ui.label("Relatedness across families").classes(
                "text-lg font-semibold text-blue-700"
            )
            ui.space()
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")
        ui.label(
            f"{len(graph.edges)} related pair(s) between {len(graph.nodes)} samples of"
            " different families, touching the individuals shown. Node size: number"
            " of relatives in other families; colour: family"
            + ("; edge colour: observed relationship." if show_detail else ".")
            + " Drag to move, scroll to zoom."
        ).classes("text-xs text-gray-500")
        if hubs:
            ui.label(
                "Hubs (red outline): "
                + ", ".join(f"{h.sample} ({h.degree})" for h in hubs[:10])
                + ". A sample related to many families at third degree is more"
                " often contaminated, or from an endogamous population, than"
                " related to them all."
            ).classes("text-xs text-purple-800")
        with ui.row().classes("w-full gap-4 no-wrap items-start"):
            ui.echart(
                cross_family_graph_option(
                    graph,
                    show_detail=show_detail,
                    highlight=keys,
                    hub_minimum=HUB_MINIMUM,
                )
            ).classes("min-w-0").style("flex: 3; height: 70vh")
            with ui.column().classes("min-w-0").style("flex: 2"):
                _cross_pair_table(graph.edges, cohort_name, show_detail)
    dialog.open()


def _render_plots(
    rows: List[Dict[str, Any]],
    by_key: Dict[Key, List[FlagRow]],
    pairs: List[PairRow],
    show_detail: bool,
) -> None:
    def flagged(row: Dict[str, Any], categories: Any) -> bool:
        key = (row.get("family_id") or "", row.get("sample") or "")
        return any(f.info.category in categories for f in by_key.get(key, []))

    sex_points = []
    for r in rows:
        point = sample_point(r)
        point["flagged"] = flagged(r, SEX_CATEGORIES)
        if not show_detail:
            point["inferred_sex"] = point["karyotype"] = None
        sex_points.append(point)
    _subtitle("Sex chromosome copy numbers")
    ui.echart(
        sex_scatter_option(sex_points, title="", show_detail=show_detail)
    ).classes("w-full h-72")
    ui.label(
        "Coloured by pedigree sex; diamonds: sex or karyotype flags. Shaded: within"
        f" ±{THRESHOLDS.get('sex_tolerance', 0.25):g} of a karyotype."
    ).classes("text-xs text-gray-500")

    if pairs:
        _subtitle("Relatedness")
        ui.echart(
            relatedness_scatter_option(
                pairs if show_detail else [p for p in pairs if p.declared],
                title="",
                show_detail=show_detail,
                label_mismatches=sum(p.status == "mismatch" for p in pairs) <= 25,
            )
        ).classes("w-full h-80")
        ui.label(
            "Within-family pairs of the shown individuals, and their relatives in"
            " other families (the cohort keeps only cross-family pairs observed"
            " related). Red diamonds: mismatches with the pedigree."
        ).classes("text-xs text-gray-500")

    charr_max = THRESHOLDS.get("charr_max", 0.05)
    charr_flag = [r for r in rows if flagged(r, QUALITY_CATEGORIES)]
    charr_rest = [r for r in rows if not flagged(r, QUALITY_CATEGORIES)]
    _subtitle("Contamination")
    ui.echart(
        metric_strip_option(
            [
                {
                    "name": r.get("sample") or "",
                    "value": to_float(r.get("contamination_charr")),
                    "flagged": True,
                }
                for r in charr_flag
            ],
            background=[
                v
                for v in (to_float(r.get("contamination_charr")) for r in charr_rest)
                if v is not None
            ],
            title="",
            axis_name="CHARR",
            lines=[(charr_max, f"flag > {charr_max:g}", "#c62828")],
        )
    ).classes("w-full h-48")

    fence = mendelian_fence(rows)
    mendel = [
        (r, to_float(r.get("trio_mendelian_error_rate")))
        for r in rows
        if to_float(r.get("trio_mendelian_error_rate")) is not None
    ]
    if mendel:
        out = [(r, v) for r, v in mendel if fence is not None and v > fence]
        rest = [v for r, v in mendel if fence is None or v <= fence]
        _subtitle("Trio Mendelian error rate")
        ui.echart(
            metric_strip_option(
                [
                    {"name": r.get("sample") or "", "value": v, "flagged": True}
                    for r, v in out
                ],
                background=rest,
                title="",
                axis_name="errors / informative sites (children)",
                lines=[(fence, "Q3 + 1.5 IQR", "#455a64")] if fence is not None else [],
            )
        ).classes("w-full h-48")
        ui.label(
            "Not flagged by the pipeline: the diamonds are outliers of the shown"
            " children, worth a look for a subtle swap or contamination."
        ).classes("text-xs text-gray-500")


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------


def _download_menu(
    store: Any,
    cohort: Any,
    state: Dict[str, Any],
    filtered_state: Dict[str, Any],
    show_detail: bool,
) -> None:
    """Flags, samples and pairs of the shown individuals, apart or zipped."""

    def keys() -> Set[Key]:
        return _keys(filtered_state)

    def flags() -> List[FlagRow]:
        qc: Optional[QcData] = state["qc"]
        if qc is None:
            return []
        shown = keys()
        return [f for f in qc.flags if any(k in shown for k in f.keys)]

    def build_flags() -> str:
        rows, columns = flags_export(flags(), load_acks(store.data_dir), show_detail)
        return rows_to_tsv(rows, columns)

    def build_samples() -> str:
        qc: Optional[QcData] = state["qc"]
        if qc is None:
            return ""
        people = {
            (p["FID"], p["IID"]): p
            for p in people_from_individuals(filtered_state.get("individuals", []))
        }
        rows, columns = samples_export(
            qc.sample_rows(), people, flags(), load_acks(store.data_dir), show_detail
        )
        return rows_to_tsv(rows, columns)

    def build_pairs() -> str:
        qc: Optional[QcData] = state["qc"]
        if qc is None:
            return ""
        shown = keys()
        pairs = [
            r
            for r in qc.pair_rows()
            if (r.get("family_id_a") or "", r.get("sample_a")) in shown
            or (r.get("family_id_b") or r.get("family_id_a") or "", r.get("sample_b"))
            in shown
        ]
        rows, columns = pairs_export(
            pairs, flags(), load_acks(store.data_dir), show_detail
        )
        return rows_to_tsv(rows, columns)

    def filename(table: Optional[str], ext: str = "tsv") -> Callable[[], str]:
        return lambda: export_filename(
            cohort.name, "sampleqc", table, tag=state["label"], ext=ext
        )

    download_button(
        [
            Download("flags", "Flags", filename("flags"), build_flags),
            Download(
                "samples", "Sex and sample quality", filename("samples"), build_samples
            ),
            Download("pairs", "Pairs (relatedness)", filename("pairs"), build_pairs),
        ],
        mark="download-sampleqc",
        zip_filename=filename(None, "zip"),
    )


# ---------------------------------------------------------------------------
# Tab
# ---------------------------------------------------------------------------


def render_cohort_sample_qc_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort-wide Sample QC sub-tab.

    Args:
        store: DataStore instance
        cohort: Cohort instance
        filtered_state: Dict with 'individuals' key holding the filtered rows
        refreshers: List to append the refresh callback to
    """
    qc_dir = get_cohort_qc_dir(cohort)
    labels = find_qc_labels(qc_dir, cohort.name)
    if not labels:
        ui.label(f"No cohort-level sample QC results in: {qc_dir}").classes(
            "text-gray-500 italic"
        )
        return

    show_detail = can_write()
    state: Dict[str, Any] = {
        "label": labels[0],
        "qc": None,
        "loading": True,
        "include_acked": False,
    }

    async def load() -> None:
        state["loading"] = True
        label = state["label"]
        try:
            qc = await asyncio.to_thread(load_cohort_qc, cohort, label)
        except Exception:
            logger.warning("Failed to read cohort QC of %s", cohort.name, exc_info=True)
            qc = None
        if label != state["label"]:
            return  # the label changed while reading
        state["qc"] = qc
        state["loading"] = False
        render_body.refresh()

    with ui.card().classes("w-full min-w-0"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("Sample QC").classes("text-lg font-semibold text-blue-700")
            if len(labels) == 1:
                ui.label(f"somalier · {labels[0]}").classes("text-xs text-gray-500")
            ui.space()

            def on_include_acked(e: Any) -> None:
                state["include_acked"] = bool(e.value)
                render_body.refresh()

            ui.switch(
                "Include acknowledged", value=False, on_change=on_include_acked
            ).props("dense").classes("text-xs").tooltip(
                "Count the flags a curator reviewed and acknowledged"
            ).mark("cohort-qc-include-acked")
            _download_menu(store, cohort, state, filtered_state, show_detail)
            report = cohort_report_url(cohort, labels[0], get_static_prefix())
            if report:
                ui.link("somalier report ↗", report, new_tab=True).classes(
                    "text-xs"
                ).tooltip("somalier's own interactive relate report, every pair")

        if len(labels) > 1:

            def on_label(e: Any) -> None:
                state.update(label=e.value, qc=None, loading=True)
                render_body.refresh()
                ui.timer(0.05, load, once=True)

            ui.select(
                options=labels, value=labels[0], label="Sites", on_change=on_label
            ).props("outlined dense").classes("w-full mb-2")

        @ui.refreshable
        def render_body() -> None:
            qc: Optional[QcData] = state["qc"]
            if state["loading"]:
                with ui.row().classes("items-center gap-2"):
                    ui.spinner(size="sm")
                    ui.label("Reading the cohort's sample QC tables…").classes(
                        "text-xs text-gray-500"
                    )
                return
            if qc is None or qc.samples is None:
                ui.label(f"No cohort samples table for {state['label']}").classes(
                    "text-gray-500 italic"
                )
                return

            keys = _keys(filtered_state)
            rows = [
                r
                for r in qc.sample_rows()
                if ((r.get("family_id") or ""), (r.get("sample") or "")) in keys
            ]
            checked = {
                ((r.get("family_id") or ""), (r.get("sample") or "")) for r in rows
            }
            acks = load_acks(store.data_dir)
            active, acked = split_acknowledged(qc.flags, acks)
            by_key = flags_by_key(qc.flags if state["include_acked"] else active)
            acked_shown = sum(1 for f in acked if any(k in checked for k in f.keys))
            all_pairs = pair_rows(qc.pairs)
            pairs = [
                p
                for p in all_pairs
                if (p.family_a, p.sample_a) in keys or (p.family_b, p.sample_b) in keys
            ]
            cross = [p for p in pairs if p.cross_family]
            graph = cross_family_graph(cross)
            cross_counts = {k: n.degree for k, n in graph.nodes.items()}

            note = f"{len(checked)} of {len(keys)} shown individuals checked."
            if len(checked) < len(keys):
                note += " The others have no somalier extract."
            if acked_shown:
                note += f" {acked_shown} acknowledged flag(s)" + (
                    " counted." if state["include_acked"] else " not counted."
                )
            ui.label(note).classes("text-xs text-gray-500")

            _render_counts(by_key, checked)

            _subtitle("Flagged individuals")
            _render_flagged(
                by_key, checked, cross_counts, cohort.name, show_detail, acks
            )

            _subtitle("Relatives across families")
            hubs = [
                n
                for k, n in graph.nodes.items()
                if k in keys and n.degree >= HUB_MINIMUM
            ]
            shown_with = {k for k in cross_counts if k in keys}
            text = (
                f"{len(cross)} related pair(s) across families, involving"
                f" {len(shown_with)} shown individual(s)."
            )
            if hubs:
                hubs.sort(key=lambda n: -n.degree)
                text += " Hubs: " + ", ".join(
                    f"{h.sample} ({h.degree})" for h in hubs[:5]
                )
            ui.label(text).classes("text-xs text-gray-600")
            if cross:
                ui.button(
                    "Cross-family network…",
                    icon="hub",
                    on_click=lambda: _open_cross_family_dialog(
                        cross, keys, cohort.name, show_detail
                    ),
                ).props("flat dense no-caps color=primary")

            _render_plots(rows, by_key, pairs, show_detail)

        render_body()
        refreshers.append(render_body.refresh)
        ui.timer(0.05, load, once=True)
