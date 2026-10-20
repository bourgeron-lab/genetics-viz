"""Sample QC tab for the family page: the ghfc-ngs ``sampleqc`` results.

Four blocks, each a table beside the plot that reads it:

1. the issues: the family's flags, plus the cohort-level ones about its
   members (relatives in other families, sample swap candidates);
2. relationships: the kinship matrix (observed against pedigree) and the
   KING / IBS0 plane, then every pair;
3. sex and karyotype: X/Y copy numbers and chrX heterozygosity, over the
   rest of the cohort;
4. sample quality: depth, sites, contamination and Mendelian errors.

The family's own files are a few KB and read at once. The cohort's add the
cross-family flags and the plots' background; they are read off-thread, so
the tab draws first and fills in.

Wording follows the user's role, as on the ROH tab: viewers get neutral
labels, no evidence, and only the pairs the pedigree declares related.
"""

import asyncio
import logging
from dataclasses import replace
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

from nicegui import ui

from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.roh import to_float
from genetics_viz.utils.sample_qc import (
    HUB_MINIMUM,
    QUALITY_CATEGORIES,
    SEX_CATEGORIES,
    THRESHOLDS,
    FlagRow,
    PairRow,
    QcData,
    cross_family_graph,
    family_flags_with_cohort,
    flag_badges,
    flag_title,
    flag_tooltip,
    flags_by_sample,
    get_family_qc_dir,
    is_typical_karyotype,
    karyotype_color,
    load_cohort_qc,
    load_family_qc,
    mendelian_fence,
    pair_rows,
    probe_cohort_qc_data,
    probe_qc_data,
    relationship_info,
    sample_point,
    severity_info,
    status_info,
)
from genetics_viz.utils.sample_qc_plots import (
    kinship_matrix_option,
    metric_strip_option,
    relatedness_scatter_option,
    sex_color,
    sex_scatter_option,
    x_het_option,
)

logger = logging.getLogger(__name__)

__all__ = ["family_qc_summary", "probe_qc_data", "render_sample_qc_tab"]

_NUM = {"cellType": "number", "sorting": "numerical", "sortable": True}


def _section(title: str, subtitle: str = "") -> None:
    ui.label(title).classes("text-base font-semibold text-blue-700 mt-4")
    if subtitle:
        ui.label(subtitle).classes("text-xs text-gray-500 mb-1")


def _table_box() -> ui.element:
    """A block box the table can scroll inside (see cohort_roh_tab)."""
    return ui.element("div").classes("w-full min-w-0")


# ---------------------------------------------------------------------------
# 1. Issues
# ---------------------------------------------------------------------------


def _render_issues(
    flags: List[FlagRow],
    family_id: str,
    cohort_name: str,
    show_detail: bool,
    cohort_note: str,
) -> None:
    if not flags:
        with ui.row().classes("items-center gap-2"):
            ui.icon("check_circle", color="green")
            ui.label("No sample QC issue among the selected members.").classes(
                "text-green-800"
            )
    else:
        rows = []
        for i, f in enumerate(flags):
            sev = severity_info(f.severity)
            other = f.other_family(family_id)
            rows.append(
                {
                    "key": str(i),
                    "severity": sev["label"],
                    "severity_color": sev["color"],
                    "severity_rank": sev["rank"],
                    "category": f.info.category_label,
                    "category_color": f.info.category_color,
                    "flag": flag_title(f, show_detail),
                    "flag_badges": [
                        {
                            "label": flag_title(f, show_detail),
                            "color": sev["color"],
                            "tooltip": flag_tooltip(f, show_detail),
                        }
                    ],
                    "samples": " – ".join(f.samples),
                    "other_family": other or "",
                    "evidence": f.evidence,
                }
            )
        columns: List[Dict[str, Any]] = [
            {
                "id": "severity",
                "header": "Severity",
                "cellType": "badge",
                "colorField": "severity_color",
                "sortable": True,
                "sortField": "severity_rank",
                "sorting": "numerical",
            },
            {
                "id": "category",
                "header": "Category",
                "cellType": "badge",
                "colorField": "category_color",
                "sortable": True,
            },
            {"id": "flag", "header": "Issue", "cellType": "badge_list"},
            {"id": "samples", "header": "Samples", "sortable": True},
        ]
        if any(r["other_family"] for r in rows) and cohort_name:
            columns.append(
                {
                    "id": "other_family",
                    "header": "Other family",
                    "cellType": "link",
                    "href": f"/cohort/{cohort_name}/family/{{other_family}}",
                    "sortable": True,
                }
            )
        if show_detail:
            columns.append({"id": "evidence", "header": "Evidence"})
        with _table_box():
            DataTable(
                columns=columns,
                rows=rows,
                row_key="key",
                pagination={"rowsPerPage": 10},
            )
        ui.label("Hover an issue for what it means and what to check.").classes(
            "text-xs text-gray-500"
        )
    if cohort_note:
        ui.label(cohort_note).classes("text-xs text-gray-500 italic")


# ---------------------------------------------------------------------------
# 2. Relationships
# ---------------------------------------------------------------------------


def _render_relationships(
    pairs: List[PairRow],
    members: List[str],
    cross: List[PairRow],
    show_detail: bool,
    relative_links: Optional[Dict[Tuple[str, str], int]] = None,
) -> None:
    within = [p for p in pairs if p.sample_a in members and p.sample_b in members]
    if not within and not cross:
        ui.label("No pair among the selected members.").classes("text-gray-500 italic")
        return
    hidden = 0
    if not show_detail:
        # Viewers get the pedigree-declared pairs only, as on the ROH tab; the
        # cross-family relatives are named, neutrally, among the issues.
        visible = [p for p in within if p.declared]
        hidden = len(within) - len(visible) + len(cross)
        within, cross = visible, []
        if not within:
            ui.label(
                "No pedigree-declared relative pairs among the selected members;"
                f" {hidden} other pair(s) are not shown."
            ).classes("text-gray-500 italic")
            return

    with ui.card().classes("w-full"):
        with ui.row().classes("w-full gap-4 no-wrap"):
            if len(members) > 1:
                ui.echart(
                    kinship_matrix_option(members, within, show_detail=show_detail)
                ).classes("min-w-0").style(
                    f"flex: 1; height: {max(300, 60 + 42 * len(members))}px"
                )
            ui.echart(
                relatedness_scatter_option(
                    within + cross,
                    title="KING kinship against IBS0/IBS2",
                    show_detail=show_detail,
                )
            ).classes("min-w-0 h-80").style("flex: 1")

    rows = []
    for p in within + cross:
        status = status_info(p.status)
        rows.append(
            {
                "key": p.key,
                "pair": f"{p.sample_a} – {p.sample_b}",
                "other_family": p.family_b if p.cross_family else "",
                "expected": relationship_info(p.expected)["label"],
                "expected_color": relationship_info(p.expected)["color"],
                "observed": relationship_info(p.observed)["label"],
                "observed_color": relationship_info(p.observed)["color"],
                "kinship": p.kinship,
                "expected_kinship": p.expected_kinship,
                "ratio": p.ratio,
                "n_sites": p.n_sites,
                "po_rate": p.po_rate,
                "status": status["label"],
                "status_color": status["color"],
                "mismatch_rank": 1 if p.status == "mismatch" else 0,
                "relative_links": (
                    (relative_links or {}).get((p.family_b, p.sample_b))
                    if p.cross_family
                    else None
                ),
            }
        )
    # Mismatches first, matching the sort indicator.
    rows.sort(key=lambda r: -r["mismatch_rank"])
    columns: List[Dict[str, Any]] = [
        {"id": "pair", "header": "Pair", "sortable": True},
        {
            "id": "expected",
            "header": "Pedigree",
            "cellType": "badge",
            "colorField": "expected_color",
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
        {"id": "expected_kinship", "header": "Pedigree φ", **_NUM},
        {"id": "ratio", "header": "IBS0/IBS2", **_NUM},
        {"id": "n_sites", "header": "Sites", **_NUM},
        {"id": "po_rate", "header": "PO opposite hom", **_NUM},
        {
            "id": "status",
            "header": "Check",
            "cellType": "badge",
            "colorField": "status_color",
            "sortable": True,
            "sortField": "mismatch_rank",
            "sorting": "numerical",
        },
    ]
    hubs = False
    if any(r["other_family"] for r in rows):
        columns.insert(1, {"id": "other_family", "header": "Other family"})
        columns.insert(
            2,
            {"id": "relative_links", "header": "Relative's cross-family links", **_NUM},
        )
        hubs = any((r["relative_links"] or 0) >= HUB_MINIMUM for r in rows)
    with _table_box():
        DataTable(
            columns=columns,
            rows=rows,
            row_key="key",
            pagination={"rowsPerPage": 10},
            initial_sorting=[{"id": "status", "desc": True}],
        )
    note = (
        "KING kinship: about 0.25 for parent–offspring and full sibs, 0.125 for"
        " second degree, 0.0625 for first cousins. Parent–offspring share an allele"
        " at every site, so their IBS0/IBS2 is near 0; full sibs sit to its right."
        " PO opposite hom: 0/0 against 1/1 in the family VCF, for declared"
        " parent–child pairs."
    )
    if hubs:
        note += (
            " A relative linked to many families is a hub: contamination or"
            " population structure is likelier than kinship with each of them."
        )
    if hidden:
        note += f" {hidden} pair(s) not declared related in the pedigree are not shown."
    ui.label(note).classes("text-xs text-gray-500")


# ---------------------------------------------------------------------------
# 3. Sex and karyotype
# ---------------------------------------------------------------------------


def _render_sex(
    rows: List[Dict[str, Any]],
    by_sample: Dict[str, List[FlagRow]],
    background: List[Tuple[float, float]],
    background_het: List[Tuple[float, float]],
    show_detail: bool,
) -> None:
    table = []
    points = []
    for r in rows:
        sample = r.get("sample") or ""
        flags = [
            f for f in by_sample.get(sample, []) if f.info.category in SEX_CATEGORIES
        ]
        karyotype = r.get("karyotype") or "unknown"
        if not show_detail and not is_typical_karyotype(karyotype):
            karyotype_shown = "Atypical" if karyotype != "unknown" else "unknown"
        else:
            karyotype_shown = karyotype
        table.append(
            {
                "sample": sample,
                "pedigree_sex": r.get("pedigree_sex") or "unknown",
                "pedigree_sex_color": sex_color(r.get("pedigree_sex")),
                "inferred_sex": r.get("inferred_sex") or "unknown",
                "inferred_sex_color": sex_color(r.get("inferred_sex")),
                "karyotype": karyotype_shown,
                "karyotype_color": karyotype_color(karyotype),
                "x_copies": to_float(r.get("x_copies")),
                "y_copies": to_float(r.get("y_copies")),
                "x_het": to_float(r.get("x_het_rate")),
                "x_sites": r.get("x_sites"),
                "evidence": r.get("sex_evidence") or "—",
                "flags": "",
                "flags_badges": flag_badges(flags, show_detail),
            }
        )
        point = sample_point(r)
        point["flagged"] = bool(flags)
        if not show_detail:
            point["inferred_sex"] = point["karyotype"] = None
        points.append(point)

    columns: List[Dict[str, Any]] = [
        {"id": "sample", "header": "Sample", "sortable": True},
        {
            "id": "pedigree_sex",
            "header": "Pedigree sex",
            "cellType": "badge",
            "colorField": "pedigree_sex_color",
            "sortable": True,
        },
        {
            "id": "inferred_sex",
            "header": "Inferred sex",
            "cellType": "badge",
            "colorField": "inferred_sex_color",
            "sortable": True,
        },
        {
            "id": "karyotype",
            "header": "Karyotype",
            "cellType": "badge",
            "colorField": "karyotype_color",
            "sortable": True,
        },
        {"id": "x_copies", "header": "X copies", **_NUM},
        {"id": "y_copies", "header": "Y copies", **_NUM},
        {"id": "x_het", "header": "X het", **_NUM},
        {"id": "x_sites", "header": "X sites", **_NUM},
        {"id": "evidence", "header": "Copy numbers from", "sortable": True},
        {"id": "flags", "header": "Flags", "cellType": "badge_list"},
    ]
    with _table_box():
        DataTable(
            columns=columns,
            rows=table,
            row_key="sample",
            pagination={"rowsPerPage": 10},
        )
    with ui.card().classes("w-full"):
        with ui.row().classes("w-full gap-4 no-wrap"):
            ui.echart(
                sex_scatter_option(
                    points, background=background, show_detail=show_detail
                )
            ).classes("min-w-0 h-80").style("flex: 1")
            ui.echart(
                x_het_option(points, background=background_het, show_detail=show_detail)
            ).classes("min-w-0 h-80").style("flex: 1")
    ui.label(
        "Copy numbers from the coverage bigWigs when the cnv_coverage step ran"
        ' ("coverage"), otherwise from somalier\'s sites. Shaded boxes: within'
        f" ±{THRESHOLDS.get('sex_tolerance', 0.25):g} of a karyotype; outside every"
        " box the call is ambiguous (mosaicism or an artefact). Points are coloured"
        " by pedigree sex; grey: the rest of the cohort."
    ).classes("text-xs text-gray-500")
    if show_detail and any(
        not is_typical_karyotype(r.get("karyotype")) and r.get("karyotype") != "unknown"
        for r in rows
    ):
        ui.label(
            "To confirm a sex chromosome anomaly, open an SV on chrX or chrY of this"
            " sample and look at its Copy Number View, then request a karyotype."
            " Note the CNV control panel has no group for XXY or X0: the copy number"
            " there is the sample's own normalised coverage."
        ).classes("text-xs text-purple-800")


# ---------------------------------------------------------------------------
# 4. Sample quality
# ---------------------------------------------------------------------------


def _render_quality(
    rows: List[Dict[str, Any]],
    by_sample: Dict[str, List[FlagRow]],
    cohort_rows: List[Dict[str, Any]],
    member_set: Set[str],
    show_detail: bool,
) -> None:
    table = []
    for r in rows:
        sample = r.get("sample") or ""
        flags = [
            f
            for f in by_sample.get(sample, [])
            if f.info.category in QUALITY_CATEGORIES
        ]
        table.append(
            {
                "sample": sample,
                "depth": to_float(r.get("depth_mean")),
                "sites": to_float(r.get("n_sites_called")),
                "het_hom": to_float(r.get("het_hom_alt_ratio")),
                "charr": to_float(r.get("contamination_charr")),
                "middling": to_float(r.get("p_middling_ab")),
                "trio_sites": to_float(r.get("trio_informative_sites")),
                "mendel": to_float(r.get("trio_mendelian_error_rate")),
                "flags": "",
                "flags_badges": flag_badges(flags, show_detail),
            }
        )
    with _table_box():
        DataTable(
            columns=[
                {"id": "sample", "header": "Sample", "sortable": True},
                {"id": "depth", "header": "Depth", **_NUM},
                {"id": "sites", "header": "Sites called", **_NUM},
                {"id": "het_hom", "header": "Het / hom alt", **_NUM},
                {"id": "charr", "header": "CHARR", **_NUM},
                {"id": "middling", "header": "Middling AB", **_NUM},
                {"id": "trio_sites", "header": "Trio sites", **_NUM},
                {"id": "mendel", "header": "Mendelian error rate", **_NUM},
                {"id": "flags", "header": "Flags", "cellType": "badge_list"},
            ],
            rows=table,
            row_key="sample",
            pagination={"rowsPerPage": 10},
        )

    others = [r for r in cohort_rows if (r.get("sample") or "") not in member_set]
    charr_max = THRESHOLDS.get("charr_max", 0.05)
    charr_points = [
        {
            "name": t["sample"],
            "value": t["charr"],
            "flagged": (t["charr"] or 0) > charr_max,
        }
        for t in table
        if t["charr"] is not None
    ]
    fence = mendelian_fence(cohort_rows or rows)
    mendel_points = [
        {
            "name": t["sample"],
            "value": t["mendel"],
            "flagged": fence is not None and (t["mendel"] or 0) > fence,
        }
        for t in table
        if t["mendel"] is not None
    ]
    with ui.card().classes("w-full"):
        with ui.row().classes("w-full gap-4 no-wrap"):
            ui.echart(
                metric_strip_option(
                    charr_points,
                    background=[
                        v
                        for v in (
                            to_float(r.get("contamination_charr")) for r in others
                        )
                        if v is not None
                    ],
                    title="Contamination (CHARR)",
                    axis_name="CHARR",
                    lines=[(charr_max, f"flag > {charr_max:g}", "#c62828")],
                )
            ).classes("min-w-0 h-56").style("flex: 1")
            if mendel_points:
                lines = (
                    [(fence, "cohort Q3 + 1.5 IQR", "#455a64")]
                    if fence is not None
                    else []
                )
                ui.echart(
                    metric_strip_option(
                        mendel_points,
                        background=[
                            v
                            for v in (
                                to_float(r.get("trio_mendelian_error_rate"))
                                for r in others
                            )
                            if v is not None
                        ],
                        title="Trio Mendelian error rate (children)",
                        axis_name="errors / informative sites",
                        lines=lines,
                    )
                ).classes("min-w-0 h-56").style("flex: 1")
    ui.label(
        "CHARR: somalier's contamination estimate from reference reads at"
        " homozygous-alternate sites. Mendelian error rate: impossible child"
        " genotypes given the parents', in the family VCF (GQ ≥ 20, DP ≥ 10). The"
        " pipeline does not flag a high rate by itself; the line marks the cohort's"
        " outliers. Grey: the rest of the cohort."
    ).classes("text-xs text-gray-500")


# ---------------------------------------------------------------------------
# Tab
# ---------------------------------------------------------------------------


def render_sample_qc_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the Sample QC tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name ("" on the standalone family page)
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    family_qc = load_family_qc(store.data_dir, family_id)
    if family_qc is None:
        ui.label(
            f"No sample QC results found in: {get_family_qc_dir(store.data_dir, family_id)}"
        ).classes("text-gray-500 italic")
        return

    show_detail = can_write()
    cohort = store.cohorts.get(cohort_name) if cohort_name else None
    has_cohort = cohort is not None and probe_cohort_qc_data(cohort)
    state: Dict[str, Any] = {
        "qc": family_qc,
        "cohort_qc": None,
        "cohort_loading": has_cohort,
    }

    async def load_cohort() -> None:
        label = state["qc"].label
        try:
            cohort_qc = await asyncio.to_thread(load_cohort_qc, cohort, label)
        except Exception:
            logger.warning("Failed to read cohort QC of %s", cohort_name, exc_info=True)
            cohort_qc = None
        if label != state["qc"].label:
            return  # the label changed while reading
        state["cohort_qc"] = cohort_qc
        state["cohort_loading"] = False
        render_content.refresh()

    with ui.row().classes("items-center gap-2 mb-2"):
        ui.label("Sample QC").classes("text-lg font-semibold text-blue-700")
        provenance = ui.label(f"somalier · {family_qc.label}").classes(
            "text-xs text-gray-500"
        )

    if len(family_qc.labels) > 1:

        def on_label(e: Any) -> None:
            state["qc"] = (
                load_family_qc(store.data_dir, family_id, e.value) or family_qc
            )
            provenance.text = f"somalier · {state['qc'].label}"
            state["cohort_qc"] = None
            state["cohort_loading"] = has_cohort
            render_content.refresh()
            if has_cohort:
                ui.timer(0.05, load_cohort, once=True)

        ui.select(
            options=family_qc.labels,
            value=family_qc.label,
            label="Sites",
            on_change=on_label,
        ).props("outlined dense").classes("w-64 mb-2")

    @ui.refreshable
    def render_content() -> None:
        qc: QcData = state["qc"]
        cohort_qc: Optional[QcData] = state["cohort_qc"]
        members = [m for m in selected_members["value"]]
        member_set = set(members)
        if not members:
            ui.label("No members selected").classes("text-gray-500 italic")
            return
        sample_rows = [r for r in qc.sample_rows() if r.get("sample") in member_set]
        qc_members = [r["sample"] for r in sample_rows]
        not_checked = [m for m in members if m not in set(qc_members)]

        flags = [
            f
            for f in family_flags_with_cohort(qc, cohort_qc, family_id)
            if any(s in member_set for s in f.samples)
        ]
        by_sample = flags_by_sample(flags)

        if state["cohort_loading"]:
            cohort_note = "Reading the cohort's tables for relatives in other families…"
        elif cohort is None:
            cohort_note = (
                "Relatives in other families are checked at the cohort level;"
                " open the family from its cohort to see them."
            )
        elif cohort_qc is None:
            cohort_note = "No cohort-level sample QC: relatives in other families are not checked."
        else:
            cohort_note = ""
        if not_checked:
            cohort_note = (
                f"Not checked (no somalier extract): {', '.join(not_checked)}. "
                + cohort_note
            )

        _section("Issues", "Flags raised by the pipeline, most severe first.")
        _render_issues(flags, family_id, cohort_name, show_detail, cohort_note)

        cross: List[PairRow] = []
        if cohort_qc is not None:
            cross = [
                p
                for p in pair_rows(cohort_qc.pairs)
                if p.cross_family
                and (
                    (p.family_a == family_id and p.sample_a in member_set)
                    or (p.family_b == family_id and p.sample_b in member_set)
                )
            ]
            # This family's member first in each cross pair.
            cross = [
                p
                if p.family_a == family_id
                else replace(
                    p,
                    family_a=p.family_b,
                    sample_a=p.sample_b,
                    family_b=p.family_a,
                    sample_b=p.sample_a,
                )
                for p in cross
            ]
        _section(
            "Relationships",
            "Pairwise relatedness from somalier's sites in the CRAMs, checked"
            " against the pedigree"
            + (" and, below, relatives in other families." if cross else "."),
        )
        links = (
            {
                key: node.degree
                for key, node in cross_family_graph(
                    pair_rows(cohort_qc.pairs)
                ).nodes.items()
            }
            if cross and cohort_qc is not None
            else {}
        )
        _render_relationships(
            pair_rows(qc.pairs), qc_members, cross, show_detail, links
        )

        cohort_rows = cohort_qc.sample_rows() if cohort_qc is not None else []
        others = [r for r in cohort_rows if r.get("family_id") != family_id]
        background = [
            (p["x"], p["y"])
            for p in map(sample_point, others)
            if p["x"] is not None and p["y"] is not None
        ]
        background_het = [
            (p["x"], p["het"])
            for p in map(sample_point, others)
            if p["x"] is not None and p["het"] is not None
        ]
        _section(
            "Sex and karyotype",
            "Sex inferred from the X and Y copy numbers, against the pedigree.",
        )
        _render_sex(sample_rows, by_sample, background, background_het, show_detail)

        _section(
            "Sample quality", "Depth, sites called, contamination and Mendelian errors."
        )
        _render_quality(sample_rows, by_sample, cohort_rows, member_set, show_detail)

    render_content()
    data_table_refreshers.append(render_content.refresh)
    if has_cohort:
        ui.timer(0.05, load_cohort, once=True)


def family_qc_summary(
    store: Any, family_id: str
) -> Tuple[Optional[QcData], Dict[str, List[FlagRow]]]:
    """The family's QC tables and its flags by sample, for the page header.

    The family files are a few KB; cohort-level flags are left to the tab.
    """
    try:
        qc = load_family_qc(store.data_dir, family_id)
    except Exception:
        logger.warning("Failed to read sample QC of %s", family_id, exc_info=True)
        return None, {}
    if qc is None:
        return None, {}
    return qc, flags_by_sample(qc.flags)
