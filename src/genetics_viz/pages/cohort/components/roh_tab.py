"""ROH tab for the family page: runs of homozygosity of each member.

Every member is shown against the unrelated reference genomes of their own
ancestry region - the group the pipeline compared them with - and only then
against the fixed published thresholds. The tab has five blocks:

1. a summary table, one row per member;
2. per reference group, a strip of the reference with the members on it, and
   the population-vs-recent plane that ``homozygosity_pattern`` reads;
3. the genome, with each member's segments by length class;
4. the genes inside long ROH;
5. pairwise kinship checked against the pedigree.

Wording follows the user's role: see ``config/roh.yaml``. Viewers get neutral
grades and only the pedigree-declared kinship pairs; curators also get the
relationship readings, the parents' kinship and every pair.
"""

import logging
from typing import Any, Callable, Dict, List, Optional

import polars as pl
from nicegui import ui

from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.export import (
    export_filename,
    frame_to_tsv,
    people_from_family,
    rows_to_tsv,
)
from genetics_viz.utils.gene_scoring import build_gene_badges, get_gene_scorer
from genetics_viz.utils.ancestry import Variant, get_region_color
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.roh import (
    EXCESS_PERCENTILE,
    F_BACKGROUND,
    GRADE_HELP,
    POOLED,
    class_color,
    find_roh_files,
    flag_badges,
    froh_export,
    get_reference_roh_dir,
    get_roh_dir,
    grade_badge,
    grade_for,
    kinship_export,
    kinship_rows,
    load_reference_froh,
    load_reference_hotspots,
    metric_value,
    ordinal,
    parents_from_froh,
    pattern_badge,
    probe_roh_data,
    quantile,
    read_roh_table,
    reference_values,
    segments_export,
    to_float,
)
from genetics_viz.utils.roh_plots import (
    genome_svg,
    pattern_scatter_option,
    reference_strip_option,
)

logger = logging.getLogger(__name__)

__all__ = ["probe_roh_data", "render_roh_tab"]

#: Strip metrics: column -> (toggle label, axis name).
_METRICS = {
    "froh_1.5mb": ("F_ROH", "F_ROH (ROH ≥ 1.5 Mb)"),
    "f_long": ("Long", "F long (≥ 10 cM, recent)"),
    F_BACKGROUND: ("Short + intermediate", "F short + intermediate (population)"),
}

#: The tool's floor on f_long for a "recent" reading (MIN_LONG_F_FOR_RECENT).
_MIN_LONG_F = 0.005

#: Genes with more than this fraction inside reference hotspots are hidden by
#: default: homozygous in many outbred people, they say little about a
#: recessive hypothesis.
_HOTSPOT_GENE_CUTOFF = 0.5


def _section(title: str, subtitle: str = "") -> None:
    ui.label(title).classes("text-base font-semibold text-blue-700 mt-4")
    if subtitle:
        ui.label(subtitle).classes("text-xs text-gray-500 mb-1")


def _group_label(group: Optional[str], n: Optional[int]) -> str:
    name = "pooled reference" if group == POOLED else f"{group} reference"
    return f"{name} (n={n})" if n else name


# ---------------------------------------------------------------------------
# 1. Summary
# ---------------------------------------------------------------------------


def _render_summary(rows: List[Dict], show_relationship: bool) -> None:
    table_rows = []
    for r in rows:
        pct = to_float(r.get("pct_froh_1.5mb"))
        table_rows.append(
            {
                "IID": r.get("IID"),
                "region": r.get("region") or "—",
                "region_color": get_region_color(r.get("region")),
                "compared": _group_label(r.get("ref_group"), r.get("ref_n")),
                "froh": to_float(r.get("froh_1.5mb")),
                "pct": ordinal(pct),
                "pct_value": pct,
                "f_short": to_float(r.get("f_short")),
                "f_intermediate": to_float(r.get("f_intermediate")),
                "f_long": to_float(r.get("f_long")),
                "nroh": r.get("nroh"),
                "largest": to_float(r.get("largest_roh_mb")),
                "pattern": r.get("homozygosity_pattern"),
                "pattern_badges": [pattern_badge(r.get("homozygosity_pattern"))],
                "grade": r.get("froh_category"),
                "grade_badges": [
                    grade_badge(r.get("froh_category"), show_relationship)
                ],
                "flags": "",
                "flags_badges": flag_badges(r),
                "parents_kinship": to_float(r.get("parents_kinship")),
                "expected_f": to_float(r.get("expected_f")),
            }
        )

    num = {"cellType": "number", "sorting": "numerical", "sortable": True}
    columns: List[Dict[str, Any]] = [
        {"id": "IID", "header": "Sample", "sortable": True},
        {
            "id": "region",
            "header": "Region",
            "cellType": "badge",
            "colorField": "region_color",
            "sortable": True,
        },
        {"id": "compared", "header": "Compared with", "sortable": True},
        {"id": "froh", "header": "F_ROH ≥1.5 Mb", **num},
        {
            "id": "pct",
            "header": "Percentile",
            "sortable": True,
            "sortField": "pct_value",
            "sorting": "numerical",
        },
        {"id": "f_short", "header": "F short", **num},
        {"id": "f_intermediate", "header": "F intermediate", **num},
        {"id": "f_long", "header": "F long", **num},
        {"id": "nroh", "header": "# ROH", **num},
        {"id": "largest", "header": "Longest (Mb)", **num},
        {"id": "pattern", "header": "Pattern", "cellType": "badge_list"},
        {"id": "grade", "header": "Grade", "cellType": "badge_list"},
        {"id": "flags", "header": "Flags", "cellType": "badge_list"},
    ]
    if show_relationship:
        columns += [
            {"id": "parents_kinship", "header": "Parents' kinship", **num},
            {"id": "expected_f", "header": "Expected F", **num},
        ]

    with ui.card().classes("w-full"):
        DataTable(
            columns=columns,
            rows=table_rows,
            row_key="IID",
            pagination={"rowsPerPage": 10},
        )
        ui.label(
            "Percentile: among unrelated reference genomes of the same ancestry region."
            " Pattern: the person against that region (hover for details). "
            + GRADE_HELP
        ).classes("text-xs text-gray-500")


# ---------------------------------------------------------------------------
# 2. Against the reference
# ---------------------------------------------------------------------------


def _render_reference_plots(
    rows: List[Dict],
    reference: Optional[pl.DataFrame],
    metric: str,
    show_relationship: bool,
) -> None:
    groups: Dict[str, List[Dict]] = {}
    for r in rows:
        groups.setdefault(r.get("ref_group") or POOLED, []).append(r)

    for group, members in sorted(groups.items()):
        ref_metric = (
            reference_values(reference, group, metric) if reference is not None else []
        )
        ref_bg = (
            reference_values(reference, group, F_BACKGROUND)
            if reference is not None
            else []
        )
        ref_long = (
            reference_values(reference, group, "f_long")
            if reference is not None
            else []
        )
        q = EXCESS_PERCENTILE / 100
        p95 = quantile(ref_metric, q)
        x_line = quantile(ref_bg, q)
        y_line = quantile(ref_long, q)
        if y_line is not None:
            y_line = max(y_line, _MIN_LONG_F)
        n = members[0].get("ref_n") or (len(ref_metric) or None)
        label = _group_label(group, n)

        strip_members = []
        scatter_points = []
        for m in members:
            color = grade_for(m.get("froh_category")).color
            value = metric_value(m, metric)
            if value is not None:
                strip_members.append({"name": m["IID"], "value": value, "color": color})
            bg, long = metric_value(m, F_BACKGROUND), to_float(m.get("f_long"))
            if bg is not None and long is not None:
                scatter_points.append(
                    {"name": m["IID"], "x": bg, "y": long, "color": color}
                )

        with ui.card().classes("w-full"):
            with ui.row().classes("w-full gap-4 no-wrap"):
                ui.echart(
                    reference_strip_option(
                        ref_metric,
                        strip_members,
                        title=f"{_METRICS[metric][1]} — {label}",
                        metric_label=_METRICS[metric][1],
                        p95=p95,
                        show_thresholds=metric == "froh_1.5mb",
                        show_relationship=show_relationship,
                    )
                ).classes("flex-1 min-w-0 h-72")
                ui.echart(
                    pattern_scatter_option(
                        list(zip(ref_bg, ref_long)),
                        scatter_points,
                        title=f"Population vs recent — {label}",
                        x_line=x_line,
                        y_line=y_line,
                    )
                ).classes("flex-1 min-w-0 h-72")


# ---------------------------------------------------------------------------
# 3. Genome
# ---------------------------------------------------------------------------


def _render_genome(
    rows: List[Dict],
    segments: Optional[pl.DataFrame],
    hotspots: Optional[pl.DataFrame],
) -> None:
    lanes = [r["IID"] for r in rows]
    if segments is None or len(segments) == 0:
        ui.label("No ROH segments for the selected members.").classes(
            "text-gray-500 italic"
        )
        return

    seg_rows = segments.filter(pl.col("IID").is_in(lanes)).to_dicts()
    by_iid: Dict[str, List] = {}
    hot_group: Dict[str, str] = {}
    for s in seg_rows:
        iid = s["IID"]
        hot_group.setdefault(iid, s.get("hotspot_group") or POOLED)
        frac = to_float(s.get("hotspot_fraction")) or 0.0
        tip = (
            f"{iid} · chr{s['chrom']}:{s['start']:,}-{s['end']:,}"
            f" · {to_float(s.get('length_mb')) or 0:.2f} Mb"
            f" / {to_float(s.get('length_cm')) or 0:.2f} cM"
            f" · {s.get('roh_class')} · {frac:.0%} in reference hotspots"
        )
        by_iid.setdefault(iid, []).append(
            (s["chrom"], int(s["start"]), int(s["end"]), s.get("roh_class"), tip)
        )

    hot_by_iid: Dict[str, List] = {}
    if hotspots is not None:
        for r in rows:
            group = hot_group.get(r["IID"]) or r.get("ref_group") or POOLED
            hot_by_iid[r["IID"]] = [
                (h["chrom"], int(h["start"]), int(h["end"]))
                for h in hotspots.filter(pl.col("group") == group).to_dicts()
            ]
    upd = {
        r["IID"]: r["upd_candidate_chrom"] for r in rows if r.get("upd_candidate_chrom")
    }

    with ui.card().classes("w-full"):
        ui.label(
            "Lanes under each chromosome, top to bottom: " + ", ".join(lanes) + "."
            " Hover a segment for its coordinates."
        ).classes("text-xs text-gray-500")
        ui.html(genome_svg(lanes, by_iid, hot_by_iid, upd), sanitize=False).classes(
            "w-full"
        )


# ---------------------------------------------------------------------------
# 4. Genes in ROH
# ---------------------------------------------------------------------------


def _render_genes(
    genes: Optional[pl.DataFrame], members: List[str], hide_hotspot: bool
) -> None:
    if genes is None:
        ui.label(
            "Gene table not produced for this family (the roh step writes it only"
            " when annotation_gencode is set)."
        ).classes("text-gray-500 italic")
        return
    rows = []
    for g in genes.filter(pl.col("IID").is_in(members)).to_dicts():
        hot = to_float(g.get("gene_hotspot_fraction")) or 0.0
        if hide_hotspot and hot > _HOTSPOT_GENE_CUTOFF:
            continue
        rows.append(
            {
                "key": f"{g['IID']}|{g.get('gene')}|{g.get('roh_start')}",
                "IID": g["IID"],
                "gene": g.get("gene"),
                "gene_id": g.get("gene_id"),
                "roh": f"chr{g['chrom']}:{g['roh_start']:,}-{g['roh_end']:,}",
                "roh_sort": f"{int(g['chrom']):02d}:{int(g['roh_start']):012d}"
                if str(g["chrom"]).isdigit()
                else str(g["chrom"]),
                "roh_length_mb": to_float(g.get("roh_length_mb")),
                "roh_class": g.get("roh_class"),
                "class_color": class_color(g.get("roh_class")),
                "in_roh": to_float(g.get("gene_fraction_in_roh")),
                "in_hotspot": hot,
            }
        )
    if not rows:
        ui.label("No genes in ROH for the selected members.").classes(
            "text-gray-500 italic"
        )
        return
    scorer = get_gene_scorer()
    for row in rows:
        row["GeneBadges"] = build_gene_badges(
            [row["gene"]] if row["gene"] else [], scorer, cap=None
        )
    # Longest ROH first. The table's sort indicator does not reorder the
    # initial rows by itself.
    rows.sort(
        key=lambda r: (-(r["roh_length_mb"] or 0), r["roh_sort"], r["gene"] or "")
    )

    num = {"cellType": "number", "sorting": "numerical", "sortable": True}
    columns = [
        {"id": "IID", "header": "Sample", "sortable": True},
        {
            "id": "gene",
            "header": "Gene",
            "sortable": True,
            "filter": {"type": "text", "placeholder": "Filter genes..."},
            "cellType": "gene_badge",
            "badgesField": "GeneBadges",
        },
        {"id": "gene_id", "header": "Gene ID", "sortable": True},
        {"id": "roh", "header": "ROH", "sortable": True, "sortField": "roh_sort"},
        {"id": "roh_length_mb", "header": "ROH (Mb)", **num},
        {
            "id": "roh_class",
            "header": "Class",
            "cellType": "badge",
            "colorField": "class_color",
            "sortable": True,
        },
        {"id": "in_roh", "header": "Fraction of gene in ROH", **num},
        {"id": "in_hotspot", "header": "Fraction in ref. hotspots", **num},
    ]
    dt_ref: Dict[str, Any] = {"dt": None}

    def on_filter(e: Dict[str, Any]) -> None:
        text = ((e.get("filters") or {}).get("gene") or "").strip().lower()
        filtered = (
            [r for r in rows if text in (r["gene"] or "").lower()] if text else rows
        )
        if dt_ref["dt"]:
            dt_ref["dt"].update_data(filtered)

    ui.label(f"{len(rows)} gene × ROH rows").classes("text-xs text-gray-500")
    dt_ref["dt"] = DataTable(
        columns=columns,
        rows=rows,
        row_key="key",
        pagination={"rowsPerPage": 25},
        on_filter=on_filter,
        initial_sorting=[{"id": "roh_length_mb", "desc": True}],
    )


# ---------------------------------------------------------------------------
# 5. Kinship
# ---------------------------------------------------------------------------


def _render_kinship(
    kinship: Optional[pl.DataFrame],
    froh: pl.DataFrame,
    members: List[str],
    show_relationship: bool,
) -> None:
    if kinship is None:
        ui.label("No kinship table (a single sequenced member).").classes(
            "text-gray-500 italic"
        )
        return
    rows = [
        r
        for r in kinship_rows(kinship, parents_from_froh(froh))
        if r["IID1"] in members and r["IID2"] in members
    ]
    hidden = 0
    if not show_relationship:
        visible = [r for r in rows if r["is_declared"]]
        hidden = len(rows) - len(visible)
        rows = visible
    for r in rows:
        r["check"] = "Mismatch" if r["mismatch"] else ("OK" if r["is_declared"] else "")
        r["check_color"] = "#c62828" if r["mismatch"] else "#2e7d32"

    if not rows:
        ui.label(
            "No pedigree-declared relative pairs among the selected members."
        ).classes("text-gray-500 italic")
    else:
        num = {"cellType": "number", "sorting": "numerical", "sortable": True}
        DataTable(
            columns=[
                {"id": "pair", "header": "Pair", "sortable": True},
                {"id": "declared", "header": "Pedigree", "sortable": True},
                {"id": "kinship", "header": "Kinship φ", **num},
                {"id": "ibs0", "header": "IBS0", **num},
                {"id": "observed", "header": "Observed", "sortable": True},
                {
                    "id": "check",
                    "header": "Check",
                    "cellType": "badge",
                    "colorField": "check_color",
                    "sortable": True,
                },
            ],
            rows=rows,
            row_key="pair",
            pagination={"rowsPerPage": 10},
        )
    note = (
        "KING-robust kinship. Parent–offspring and full sibs sit near 0.25 (1st degree,"
        " 0.177–0.354); parent–offspring have IBS0 near 0. A mismatch with the pedigree"
        " points to a sample swap or a pedigree error."
    )
    if hidden:
        note += f" {hidden} pair(s) not declared related in the pedigree are not shown."
    ui.label(note).classes("text-xs text-gray-500")


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------


def _download_menu(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    state: Dict[str, Any],
    controls: Dict[str, Any],
    show_relationship: bool,
) -> None:
    """Summary, segments, genes in ROH and kinship of the selected members."""

    def members() -> List[str]:
        return list(selected_members["value"])

    def build_summary() -> str:
        froh: pl.DataFrame = state["tables"]["froh"]
        people = {
            p["IID"]: p
            for p in people_from_family(store, cohort_name, family_id, members())
        }
        paired = [
            (r, people[r["IID"]]) for r in froh.to_dicts() if r.get("IID") in people
        ]
        rows, columns = froh_export(paired, show_relationship)
        return rows_to_tsv(rows, columns)

    def build_segments() -> str:
        frame = segments_export(state["tables"].get("roh"), family_id, members())
        return frame_to_tsv(frame) if frame is not None else ""

    def build_genes() -> str:
        genes = state["tables"].get("roh_genes")
        if (
            genes is not None
            and controls["hide_hotspot_genes"]
            and "gene_hotspot_fraction" in genes.columns
        ):
            genes = genes.filter(
                pl.col("gene_hotspot_fraction").fill_null(0.0) <= _HOTSPOT_GENE_CUTOFF
            )
        frame = segments_export(genes, family_id, members())
        return frame_to_tsv(frame) if frame is not None else ""

    def build_kinship() -> str:
        rows, columns = kinship_export(
            state["tables"].get("kinship"),
            state["tables"]["froh"],
            members(),
            show_relationship,
        )
        return rows_to_tsv(rows, columns)

    def filename(table: Optional[str], ext: str = "tsv") -> Any:
        return lambda: export_filename(
            family_id, "roh", table, tag=state["files"].file_tag, ext=ext
        )

    download_button(
        [
            Download(
                "summary", "Summary (F_ROH, grade)", filename("summary"), build_summary
            ),
            Download("segments", "ROH segments", filename("segments"), build_segments),
            Download("genes", "Genes in ROH", filename("genes"), build_genes),
            Download("kinship", "Kinship", filename("kinship"), build_kinship),
        ],
        mark="download-roh",
        zip_filename=filename(None, "zip"),
    )


# ---------------------------------------------------------------------------
# Tab
# ---------------------------------------------------------------------------


def _load(files: Any) -> Dict[str, Optional[pl.DataFrame]]:
    """Read every table of one variant. Optional tables come back ``None``."""
    out: Dict[str, Optional[pl.DataFrame]] = {}
    for kind in ("froh", "roh", "roh_genes", "kinship"):
        try:
            out[kind] = read_roh_table(files.path(kind))
        except Exception:
            logger.warning("Failed to read %s", files.path(kind), exc_info=True)
            out[kind] = None
    return out


def render_roh_tab(
    store: Any,
    family_id: str,
    cohort_name: str,
    selected_members: Dict[str, List[str]],
    data_table_refreshers: List[Callable[[], None]],
) -> None:
    """Render the ROH tab panel content.

    Args:
        store: DataStore instance
        family_id: Family ID
        cohort_name: Cohort name
        selected_members: Dict with 'value' key containing list of selected member IDs
        data_table_refreshers: List to append refresh functions to
    """
    roh_dir = get_roh_dir(store.data_dir, family_id)
    files = find_roh_files(store.data_dir, family_id)
    if files is None or files.path("froh") is None:
        ui.label(f"No ROH results found in: {roh_dir}").classes("text-gray-500 italic")
        return

    show_relationship = can_write()
    state: Dict[str, Any] = {
        "files": files,
        "tables": {},
        "reference": None,
        "hotspots": None,
    }
    controls: Dict[str, Any] = {"metric": "froh_1.5mb", "hide_hotspot_genes": True}

    def load_variant(variant: Optional[Variant] = None) -> bool:
        selected = find_roh_files(store.data_dir, family_id, variant) or files
        state["files"] = selected
        state["tables"] = _load(selected)
        state["reference"] = load_reference_froh(selected.bundle_version)
        state["hotspots"] = load_reference_hotspots(selected.bundle_version)
        froh = state["tables"].get("froh")
        return froh is not None and len(froh) > 0

    if not load_variant():
        ui.label("No individuals in the ROH table").classes("text-gray-500 italic")
        return

    with ui.row().classes("items-center gap-2 mb-2 w-full"):
        ui.label("Runs of homozygosity").classes("text-lg font-semibold text-blue-700")
        provenance = ui.label(state["files"].label).classes("text-xs text-gray-500")
        ui.space()
        _download_menu(
            store,
            family_id,
            cohort_name,
            selected_members,
            state,
            controls,
            show_relationship,
        )

    if len(files.all_variants) > 1:

        def on_variant(event: Any) -> None:
            bundle, tag = event.value.split(" / ", 1)
            load_variant((bundle, tag))
            provenance.text = state["files"].label
            render_roh_content.refresh()

        ui.select(
            options=[f"{b} / {t}" for b, t in files.all_variants],
            value=f"{files.bundle_version} / {files.filter_tag}",
            label="Bundle / filters",
            on_change=on_variant,
        ).props("outlined dense").classes("w-64 mb-2")

    @ui.refreshable
    def render_roh_content() -> None:
        tables = state["tables"]
        froh: pl.DataFrame = tables["froh"]
        reference = state["reference"]
        members = [m for m in selected_members["value"]]
        rows = [r for r in froh.to_dicts() if r.get("IID") in members]

        if reference is None:
            ref_dir = get_reference_roh_dir(state["files"].bundle_version)
            detail = (
                "set 'ancestry_reference_dir' in the application config"
                if ref_dir is None
                else f"expected at {ref_dir}"
            )
            ui.label(
                f"ROH reference unavailable for bundle {state['files'].bundle_version}"
                f" — {detail}. Showing the family's values without the reference"
                " distributions."
            ).classes("text-orange-600 text-sm mb-2")

        if not rows:
            ui.label("No members selected").classes("text-gray-500 italic")
            return

        _render_summary(rows, show_relationship)

        _section(
            "Compared with the reference",
            "Grey: unrelated reference genomes of the same ancestry region."
            f" Dotted: the region's {EXCESS_PERCENTILE:g}th percentile. Members are"
            " coloured by grade.",
        )
        ui.toggle(
            {k: v[0] for k, v in _METRICS.items()},
            value=controls["metric"],
            on_change=lambda e: (
                controls.update({"metric": e.value}),
                render_roh_content.refresh(),
            ),
        ).props("dense no-caps")
        _render_reference_plots(rows, reference, controls["metric"], show_relationship)

        _section(
            "Genome", "ROH segments by length class; shaded: reference ROH hotspots."
        )
        _render_genome(rows, tables.get("roh"), state["hotspots"])

        _section(
            "Genes in ROH",
            "Genes overlapping a ROH of at least 1 Mb (roh_min_gene_roh_mb).",
        )
        ui.checkbox(
            f"Hide genes more than {_HOTSPOT_GENE_CUTOFF:.0%} inside reference hotspots",
            value=controls["hide_hotspot_genes"],
            on_change=lambda e: (
                controls.update({"hide_hotspot_genes": e.value}),
                render_roh_content.refresh(),
            ),
        ).props("dense")
        _render_genes(tables.get("roh_genes"), members, controls["hide_hotspot_genes"])

        _section("Kinship", "Pairwise kinship checked against the pedigree.")
        _render_kinship(tables.get("kinship"), froh, members, show_relationship)

    render_roh_content()
    data_table_refreshers.append(render_roh_content.refresh)
