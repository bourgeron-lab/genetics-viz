"""Cohort-wide De novo tab: the ``denovo`` step's per-child QC and flags.

Follows the filters of the individuals table: the children shown are those
the filters leave, and every count, chart and table is theirs. From top:

1. the cohort against the literature (``references`` in ``config/denovo.yaml``);
2. strict candidates per callable Gb per child, with the published rate and the
   QC range; paternal share and transitions per child;
3. the children's QC table, linked to their family's De novo sub-tab;
4. the flags: how many candidates each marks, and what curators said of them;
5. curated de novo calls the strict set loses, or the candidates lack.

The cohort tables and the curation files are read off-thread.
"""

from __future__ import annotations

import asyncio
import logging
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from nicegui import ui

from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils import denovo, denovo_plots
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.denovo import probe_cohort_denovo
from genetics_viz.utils.export import export_filename, rows_to_tsv

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_denovo", "render_cohort_denovo_tab"]

_NUM = {"sortable": True, "sorting": "numerical"}

_CHILD_EXPORT = [
    "FID", "child", "father", "mother", "sex", "qc", "n_candidates", "n_strict",
    "n_strict_snv", "n_strict_indel", "rate_per_gb", "snv_per_gb", "indel_per_gb",
    "indel_share", "ti_fraction", "n_phased", "paternal_share", "n_parent_alt",
    "n_child_mosaic", "n_sib_shared", "n_sib_shared_strict", "callable_gb",
    "callable_status", "qc_low", "qc_high", "qc_reference", "sampleqc_flags",
]  # fmt: skip

_CANDIDATE_EXPORT = [
    "family_id", "child", "father", "mother", "child_sex", "chrom", "pos", "ref",
    "alt", "var_type", "genes", "consequences", "strict", "flags", "child_alt",
    "child_ref", "child_vaf", "father_alt", "mother_alt", "sibs_with_alt", "origin",
    "origin_status", "recurrence_families",
]  # fmt: skip


def _box() -> ui.element:
    """A block box the table can scroll inside (see cohort_roh_tab)."""
    return ui.element("div").classes("w-full min-w-0")


def _keys(filtered_state: Dict[str, Any]) -> Set[Tuple[str, str]]:
    return {
        (str(ind.get("FID")), str(ind["Sample ID"]))
        for ind in filtered_state.get("individuals", [])
        if ind.get("Sample ID")
    }


def _family_link(cohort: Any, extra: str = "") -> str:
    return f"/cohort/{cohort.name}/family/{{FID}}?tab=snvs&wcfg=denovo{extra}"


def _load(cohort: Any, data_dir: Any, label: Optional[str]) -> Dict[str, Any]:
    """Blocking: the cohort tables and the curation of their families."""
    data = denovo.load_cohort_denovo(cohort, label)
    if data is None:
        return {"data": None}
    families = (
        {str(f) for f in data.dnm["family_id"].to_list()} if len(data.dnm) else set()
    )
    curation = denovo.load_curation(data_dir, families)
    return {"data": data, "curation": curation, "families": families}


def render_cohort_denovo_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort-wide De novo tab."""
    directory = denovo.get_cohort_denovo_dir(cohort)
    labels = denovo.find_denovo_labels(directory, cohort.name)
    if not labels:
        ui.label(f"No cohort-level denovo results in: {directory}").classes(
            "text-gray-500 italic"
        )
        return
    show_curation = can_write()
    state: Dict[str, Any] = {"label": labels[0], "loaded": None, "loading": True}

    async def load() -> None:
        state["loading"] = True
        label = state["label"]
        try:
            loaded = await asyncio.to_thread(_load, cohort, store.data_dir, label)
        except Exception:
            logger.warning(
                "Failed to read the denovo tables of %s", cohort.name, exc_info=True
            )
            loaded = {"data": None}
        if label != state["label"]:
            return
        state.update(loaded=loaded, loading=False)
        render_body.refresh()

    def shown() -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """The shown children's QC rows and candidates."""
        data = state["loaded"]["data"]
        keys = _keys(filtered_state)
        children = denovo.child_rows(data, keys)
        wanted = {(c["FID"], c["child"]) for c in children}
        candidates = [
            r
            for r in data.rows()
            if (str(r.get("family_id")), str(r.get("child"))) in wanted
        ]
        return children, candidates

    with ui.card().classes("w-full min-w-0"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("De novo SNVs and indels").classes(
                "text-lg font-semibold text-blue-700"
            )
            ui.label(f"denovo step · {labels[0] if len(labels) == 1 else ''}").classes(
                "text-xs text-gray-500"
            )
            ui.space()

            def export_children() -> str:
                if not state["loaded"] or state["loaded"]["data"] is None:
                    return ""
                return rows_to_tsv(shown()[0], _CHILD_EXPORT)

            def export_candidates() -> str:
                if not state["loaded"] or state["loaded"]["data"] is None:
                    return ""
                return rows_to_tsv(shown()[1], _CANDIDATE_EXPORT)

            download_button(
                [
                    Download(
                        "children",
                        "Per-child QC of the children shown",
                        lambda: export_filename(
                            cohort.name, "denovo", "children", tag=state["label"]
                        ),
                        export_children,
                    ),
                    Download(
                        "candidates",
                        "Every candidate of the children shown, with its flags",
                        lambda: export_filename(
                            cohort.name, "denovo", "candidates", tag=state["label"]
                        ),
                        export_candidates,
                    ),
                ],
                mark="download-denovo",
            )

        if len(labels) > 1:

            def on_label(e: Any) -> None:
                state.update(label=e.value, loaded=None, loading=True)
                render_body.refresh()
                ui.timer(0.05, load, once=True)

            ui.select(
                options=labels, value=labels[0], label="Label", on_change=on_label
            ).props("outlined dense").classes("w-64 mb-2")

        @ui.refreshable
        def render_body() -> None:
            if state["loading"]:
                with ui.row().classes("items-center gap-2"):
                    ui.spinner(size="sm")
                    ui.label("Reading the cohort's de novo tables…").classes(
                        "text-xs text-gray-500"
                    )
                return
            loaded = state["loaded"] or {}
            data = loaded.get("data")
            if data is None or data.qc is None:
                ui.label(
                    f"No cohort de novo tables for label {state['label']}"
                ).classes("text-gray-500 italic")
                return
            children, candidates = shown()
            if not children:
                ui.label(
                    "No child with de novo results among the individuals shown"
                ).classes("text-gray-500 italic")
                return
            summary = denovo.cohort_summary(children)
            masks = data.meta.get("masks")
            ui.label(
                f"{len(children)} children shown, {len(candidates):,} candidates,"
                f" {summary['n_strict']:,} strict on the autosomes"
                f" ({summary['n_snv']:,} SNVs, {summary['n_indel']:,} indels)."
                + (f" Masks {masks}." if masks else "")
                + f" Median callable unmasked genome {summary['callable_gb_median'] or 0:.2f} Gb."
            ).classes("text-xs text-gray-500")

            _render_references(summary)

            with ui.row().classes("w-full gap-4 no-wrap mt-2"):
                ui.echart(denovo_plots.rate_option(children)).classes(
                    "flex-1 min-w-0"
                ).style("height: 380px").mark("denovo-rate-chart")
                ui.echart(denovo_plots.shares_option(children)).classes(
                    "flex-1 min-w-0"
                ).style("height: 380px").mark("denovo-shares-chart")

            _render_children(cohort, children)

            curation = loaded.get("curation") if show_curation else None
            flags = denovo.flag_summary(candidates, curation)
            ui.label("Flags").classes("text-base font-semibold text-blue-700 mt-4")
            with ui.row().classes("w-full gap-4 no-wrap"):
                ui.echart(denovo_plots.flags_option(flags, len(candidates))).classes(
                    "w-[34%] min-w-0 shrink-0"
                ).style(f"height: {120 + 26 * len(flags)}px").mark("denovo-flags-chart")
                with _box():
                    _render_flag_table(flags, show_curation)

            if show_curation and curation is not None:
                _render_curation(cohort, candidates, curation, loaded["families"])

        render_body()
        refreshers.append(render_body.refresh)
        ui.timer(0.05, load, once=True)


def _fmt(fmt: str, value: Any) -> str:
    if value is None or not fmt:
        return "–"
    try:
        return fmt.format(value)
    except (ValueError, TypeError):
        return str(value)


def _render_references(summary: Dict[str, Any]) -> None:
    """This cohort beside the published values, one metric per line."""
    ui.label("Against the literature").classes(
        "text-base font-semibold text-blue-700 mt-2"
    )
    with (
        ui.element("div")
        .classes("grid gap-x-4 gap-y-1 text-sm w-full items-baseline")
        .style(
            "grid-template-columns: minmax(12rem, 20rem) minmax(6rem, 8rem) minmax(10rem, 14rem) minmax(0, 1fr)"
        )
        .mark("denovo-references")
    ):
        for h in ("Metric", "Shown children", "Literature", "Source"):
            ui.label(h).classes("text-xs font-semibold text-gray-500")
        for ref in denovo.REFERENCES:
            label = ui.label(ref["label"])
            if ref.get("note"):
                label.classes("underline decoration-dotted cursor-help").tooltip(
                    ref["note"]
                )
            ui.label(_fmt(ref.get("format", ""), summary.get(ref["metric"]))).classes(
                "font-semibold"
            )
            ui.label(ref["literature"])
            ui.label(ref["source"]).classes("text-xs text-gray-600")
    ui.label(
        "Rates are per child per Gb of callable, unmasked haploid sequence (2 × μ);"
        " the masks leave out a fifth of the genome, where the true rate is higher,"
        " so the published figures are a slight upper bound here. Medians over the"
        " children shown; shares pooled."
    ).classes("text-xs text-gray-500")


def _render_children(cohort: Any, children: List[Dict[str, Any]]) -> None:
    ui.label("Children").classes("text-base font-semibold text-blue-700 mt-4")
    rows = sorted(children, key=lambda c: (c["qc"] == "ok", -(c["rate_per_gb"] or 0)))
    with _box():
        DataTable(
            columns=[
                {
                    "id": "child",
                    "header": "Child",
                    "sortable": True,
                    "filter": {"type": "text"},
                },
                {
                    "id": "FID",
                    "header": "Family",
                    "cellType": "link",
                    "href": _family_link(cohort),
                    "sortable": True,
                },
                {
                    "id": "qc",
                    "header": "QC",
                    "cellType": "badge_list",
                    "badgesField": "QcBadges",
                    "sortable": True,
                },
                {"id": "n_candidates", "header": "Candidates", **_NUM},
                {"id": "n_strict", "header": "Strict", **_NUM},
                {"id": "n_strict_snv", "header": "SNVs", **_NUM},
                {"id": "n_strict_indel", "header": "Indels", **_NUM},
                {"id": "rate_per_gb", "header": "Per Gb", **_NUM},
                {"id": "ti_fraction", "header": "Ti", **_NUM},
                {"id": "n_phased", "header": "Phased", **_NUM},
                {"id": "paternal_share", "header": "Paternal", **_NUM},
                {"id": "n_parent_alt", "header": "Parent alt", **_NUM},
                {"id": "n_child_mosaic", "header": "Mosaic?", **_NUM},
                {"id": "n_sib_shared", "header": "Sib. shared", **_NUM},
                {"id": "callable_gb", "header": "Callable Gb", **_NUM},
                {"id": "sampleqc_flags", "header": "Sample QC flags", "sortable": True},
            ],  # fmt: skip
            rows=rows,
            row_key="id",
            pagination={"rowsPerPage": 15},
        )
    ui.label(
        "Strict, SNVs, Indels, Per Gb, Ti and Paternal: strict candidates on the"
        " autosomes (Paternal: share of the phased ones). Parent alt, Mosaic? and"
        " Sib. shared: candidates carrying the flag, strict or not. QC: the"
        " pipeline's verdict on the strict rate."
    ).classes("text-xs text-gray-500")


def _render_flag_table(flags: List[Dict[str, Any]], show_curation: bool) -> None:
    columns: List[Dict[str, Any]] = [
        {"id": "flag", "header": "Flag", "cellType": "badge_list", "badgesField": "FlagBadges", "sortable": True},
        {"id": "category", "header": "Category", "sortable": True},
        {"id": "n", "header": "Candidates", **_NUM},
        {"id": "share", "header": "Share", **_NUM},
        {"id": "alone", "header": "Only flag", **_NUM},
    ]  # fmt: skip
    if show_curation:
        columns += [
            {"id": "present", "header": "Validated de novo", **_NUM},
            {"id": "absent", "header": "Validated absent", **_NUM},
            {"id": "diagnosed", "header": "Diagnostic", **_NUM},
        ]
    DataTable(columns=columns, rows=flags, row_key="id", pagination={"rowsPerPage": 20})
    if show_curation:
        ui.label(
            "Validated de novo: present without an inherited mode. A flag that marks"
            " validated de novo calls costs real ones when it hides them."
        ).classes("text-xs text-gray-500")


def _render_curation(
    cohort: Any,
    candidates: List[Dict[str, Any]],
    curation: denovo.Curation,
    families: Set[str],
) -> None:
    shown_families = {str(r.get("family_id")) for r in candidates}
    shown_keys = {(str(r.get("family_id")), str(r.get("child"))) for r in candidates}
    issues = [
        r
        for r in denovo.curation_check(candidates, curation, families & shown_families)
        if (r["FID"], r["sample"]) in shown_keys or r["problem"] == "missing"
    ]
    table = denovo.strict_vs_curation(candidates, curation)
    ui.label("Curation against the strict set").classes(
        "text-base font-semibold text-blue-700 mt-4"
    )
    yes, no = table.get("yes", {}), table.get("no", {})
    ui.label(
        f"Validated candidates: strict {yes.get('present', 0)} de novo,"
        f" {yes.get('absent', 0)} absent; flagged {no.get('present', 0)} de novo,"
        f" {no.get('absent', 0)} absent."
    ).classes("text-sm").mark("denovo-curation-summary")
    if not issues:
        ui.label(
            "No curated de novo call is hidden by the strict set or missing from"
            " the candidates."
        ).classes("text-sm text-green-800").mark("denovo-curation-ok")
        return
    for r in issues:
        r["problem_label"] = (
            "Not strict" if r["problem"] == "hidden" else "Not a candidate"
        )
    with _box():
        DataTable(
            columns=[
                {"id": "problem_label", "header": "Problem", "sortable": True},
                {"id": "sample", "header": "Sample", "sortable": True},
                {
                    "id": "FID",
                    "header": "Family",
                    "cellType": "link",
                    "href": _family_link(cohort, "&wdm=all"),
                    "sortable": True,
                },
                {
                    "id": "Variant",
                    "header": "Variant",
                    "sortable": True,
                    "sorting": "genomic",
                },
                {"id": "genes", "header": "Genes", "sortable": True},
                {"id": "validation", "header": "Validation", "sortable": True},
                {"id": "inheritance", "header": "Inheritance", "sortable": True},
                {"id": "diagnostic", "header": "Diagnostic", "sortable": True},
                {
                    "id": "flags",
                    "header": "Flags",
                    "cellType": "badge_list",
                    "badgesField": "FlagBadges",
                },
            ],  # fmt: skip
            rows=issues,
            row_key="id",
            pagination={"rowsPerPage": 10},
        )
    ui.label(
        "Not strict: a candidate validated present without an inherited mode, or with"
        " a diagnostic, that a flag keeps out of the strict set. Not a candidate: a"
        " validation saying de novo with no row in the Wombat de novo table."
    ).classes("text-xs text-gray-500")
