"""Gene and variant pages: what we have seen of a gene, or a variant, anywhere.

``/gene/<symbol>`` gathers the gene's SNVs, SVs, STR genotypes, ROH segments and
diagnostics from every cohort of every data directory, through the gene index
(``utils/gene_index.py``). ``/variant/<chr:pos:ref:alt>`` answers "have we seen
this variant?": its carriers, and every validation and diagnostic recorded for
it. Both only read; curation stays on the family pages they link to.

``/goto/<n>`` switches the user's data directory to ``/data-<n>``'s before
opening a family of it, since a family page reads the user's current one.
"""

import asyncio
from typing import Any, Dict, List, Optional, Set

from fastapi import Request
from fastapi.responses import RedirectResponse
from nicegui import app, ui

from genetics_viz.components.header import create_header
from genetics_viz.components.page_link import PageLink, safe_next
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils import gene_index
from genetics_viz.utils.auth import check_auth
from genetics_viz.utils.data import get_static_prefix_map
from genetics_viz.utils.gene_rows import (
    data_dir_infos,
    diagnostic_rows,
    gene_keys,
    load_curation,
    roh_rows,
    seen_summary,
    snv_carrier_rows,
    snv_variant_rows,
    str_genotyped,
    str_rows,
    sv_rows,
    validation_history,
)
from genetics_viz.utils.gene_scoring import gene_badge, gene_badge_html, get_gene_scorer
from genetics_viz.utils.genetrek import gene_url, get_snapshot

#: Gene page tabs: link slug -> label.
GENE_TABS = {
    "snv": "SNVs",
    "sv": "SVs",
    "str": "STRs",
    "roh": "ROH",
    "diag": "Diagnostics",
}

_FAMILY_LINK = {"cellType": "link", "href": "/goto/{_dir}?next={_next}"}
_VARIANT_LINK = {"cellType": "link", "href": "/variant/{Variant}"}


# ---------------------------------------------------------------------------
# Data directory switch
# ---------------------------------------------------------------------------


@ui.page("/goto/{dir_index}")
def goto_page(dir_index: str, next: str = "/") -> Optional[RedirectResponse]:
    """Make ``/data-<dir_index>`` the user's data directory, then open ``next``."""
    if redirect := check_auth():
        return redirect
    prefix = f"/data-{dir_index}"
    for path_str, p in get_static_prefix_map().items():
        if p == prefix:
            app.storage.user["data_dir"] = path_str
            break
    return RedirectResponse(safe_next(next) or "/")


# ---------------------------------------------------------------------------
# Shared pieces
# ---------------------------------------------------------------------------


def _col(id_: str, header: Optional[str] = None, **extra: Any) -> Dict[str, Any]:
    return {"id": id_, "header": header or id_, "sortable": True, **extra}


def _person_columns(multi_dir: bool) -> List[Dict[str, Any]]:
    cols = [_col("Data dir")] if multi_dir else []
    return cols + [
        _col("Cohort"),
        _col("FID", "Family", **_FAMILY_LINK),
        _col("Sample"),
    ]


def _curation_columns() -> List[Dict[str, Any]]:
    return [
        _col("Validation", cellType="validation"),
        _col("Diagnostic", cellType="badge_list", badgesField="DiagnosticBadges"),
    ]


def _gene_column() -> Dict[str, Any]:
    return _col("Genes", "Gene", cellType="gene_badge", badgesField="GeneBadges")


def _empty(text: str) -> None:
    ui.label(text).classes("text-gray-500 italic")


def _seen_text(seen: Dict[str, int]) -> str:
    def n(count: int, one: str, many: str) -> str:
        return f"{count} {one if count == 1 else many}"

    parts = [
        n(seen["individuals"], "individual", "individuals"),
        n(seen["families"], "family", "families"),
        n(seen["cohorts"], "cohort", "cohorts"),
    ]
    if seen["data_dirs"] > 1:
        parts.append(n(seen["data_dirs"], "data directory", "data directories"))
    return ", ".join(parts)


def _index_line(on_finished: Any) -> None:
    """The index's state; while it builds, re-run ``on_finished`` once it is done."""
    label = ui.label().classes("text-xs text-gray-500").mark("gene-index-status")

    def update() -> bool:
        st = gene_index.index_status()
        if st.building:
            label.set_text(
                f"Indexing cohort results: {st.done} of {st.total or '?'} tables "
                "checked. Results may be incomplete; they update when it is done."
            )
            return True
        when = st.checked_at.strftime("%H:%M") if st.checked_at else "never"
        text = f"{gene_index.indexed_tables()} cohort tables indexed, checked {when}."
        if st.errors:
            text += f" {len(st.errors)} could not be read."
            label.tooltip("\n".join(f"{p}: {e}" for p, e in st.errors.items()))
        label.set_text(text)
        return False

    if update():

        def poll() -> None:
            if not update():
                timer.cancel()
                on_finished()

        timer = ui.timer(2.0, poll)


# ---------------------------------------------------------------------------
# Gene page
# ---------------------------------------------------------------------------


def _load_gene(keys: Set[str], scorer: Any) -> Dict[str, Any]:
    """Every table's rows for the gene (blocking: Parquet scans, curation files)."""
    dirs = data_dir_infos()
    curation = load_curation(dirs)
    frames = {kind: gene_index.lookup_gene(kind, keys) for kind in gene_index.KINDS}
    carriers = snv_carrier_rows(frames["snv"], curation, dirs, scorer)
    return {
        "multi_dir": len(dirs) > 1,
        "snv_carriers": carriers,
        "snv_variants": snv_variant_rows(carriers),
        "sv": sv_rows(frames["sv"], curation, dirs, scorer),
        "str": str_rows(frames["str"], curation, dirs, scorer),
        "str_genotyped": str_genotyped(frames["str"]),
        "roh": roh_rows(frames["roh"], dirs),
        "diag": diagnostic_rows(dirs, dirs, scorer, gene_keys=keys),
    }


@ui.page("/gene/{symbol}")
def gene_page(symbol: str, request: Request) -> Optional[RedirectResponse]:
    """Everything seen of one gene, across cohorts and data directories."""
    if redirect := check_auth():
        return redirect
    link = PageLink.from_request(request)
    create_header()
    snapshot = get_snapshot()
    current, keys = gene_keys(symbol, snapshot)
    scorer = get_gene_scorer()

    with ui.column().classes("w-full px-6 py-6 gap-3"):
        with ui.row().classes("items-center gap-3"):
            ui.label(f"🧬 {current}").classes("text-3xl font-bold text-blue-900")
            ui.html(gene_badge_html(gene_badge(current, scorer)), sanitize=False)
            ui.link("GeneTrek", gene_url(current), new_tab=True).classes("text-sm")
        if current.upper() != symbol.strip().upper():
            ui.label(f"{symbol} is now {current}; both names are searched.").classes(
                "text-sm text-gray-600"
            )
        _gene_info(current, snapshot)

        content = ui.column().classes("w-full gap-3")

        async def load() -> None:
            content.clear()
            with content:
                ui.spinner(size="lg")
            data = await asyncio.to_thread(_load_gene, keys, scorer)
            content.clear()
            with content:
                _render_gene(data, link)

        _index_line(lambda: ui.timer(0, load, once=True))
        ui.timer(0.1, load, once=True)


def _gene_info(current: str, snapshot: Any) -> None:
    """GeneTrek's name, values, modes of inheritance and lists for the gene."""
    info = snapshot.info.get(current)
    if info is not None:
        if info.name:
            ui.label(info.name).classes("text-gray-700")
        facts = [f"{label} {value}" for label, value in info.values]
        modes = sorted({mode for _, mode in info.inheritance})
        if modes:
            facts.append("Inheritance: " + ", ".join(modes))
        if facts:
            ui.label(" · ".join(facts)).classes("text-sm text-gray-600")
    sets = [
        snapshot.geneset(name) for name in sorted(snapshot.members.get(current, ()))
    ]
    labels = sorted({gs.label for gs in sets if gs is not None})
    if labels:
        with ui.row().classes("gap-1 flex-wrap"):
            for label in labels:
                ui.label(label).classes(
                    "text-xs bg-gray-100 text-gray-700 rounded px-2 py-0.5"
                )
    if not snapshot.ok:
        ui.label(snapshot.status_text()).classes("text-xs text-gray-500")


def _render_gene(data: Dict[str, Any], link: PageLink) -> None:
    carriers = data["snv_carriers"]
    counts = {
        "snv": len(data["snv_variants"]),
        "sv": len(data["sv"]),
        "str": len(data["str"]),
        "roh": len(data["roh"]),
        "diag": len(data["diag"]),
    }
    if carriers:
        ui.label(
            f"{counts['snv']} small variant{'s' * (counts['snv'] != 1)} in "
            + _seen_text(seen_summary(carriers))
        ).classes("text-base text-gray-800").mark("gene-summary")

    linked = link.get("tab")
    first_with_rows = next((slug for slug, n in counts.items() if n), "snv")
    initial = linked if linked in GENE_TABS else first_with_rows
    multi = data["multi_dir"]

    with ui.tabs().classes("w-full").mark("gene-tabs") as tabs:
        for slug, label in GENE_TABS.items():
            ui.tab(label, label=f"{label} · {counts[slug]}")

    def on_tab(e: Any) -> None:
        slug = next((s for s, lbl in GENE_TABS.items() if lbl == e.value), None)
        link.set("tab", slug, default=first_with_rows)
        link.push()

    with ui.tab_panels(tabs, value=GENE_TABS[initial], on_change=on_tab).classes(
        "w-full"
    ):
        with ui.tab_panel(GENE_TABS["snv"]).classes("p-0 gap-4"):
            _snv_tables(data["snv_variants"], carriers, multi)
        with ui.tab_panel(GENE_TABS["sv"]).classes("p-0"):
            _sv_table(data["sv"], multi)
        with ui.tab_panel(GENE_TABS["str"]).classes("p-0"):
            _str_table(data["str"], data["str_genotyped"], multi)
        with ui.tab_panel(GENE_TABS["roh"]).classes("p-0"):
            _roh_table(data["roh"], multi)
        with ui.tab_panel(GENE_TABS["diag"]).classes("p-0"):
            _diagnostic_table(data["diag"], multi)


def _snv_tables(
    variants: List[Dict[str, Any]], carriers: List[Dict[str, Any]], multi: bool
) -> None:
    if not variants:
        _empty("No small variant of this gene in any cohort Wombat table.")
        return
    ui.label("Variants").classes("text-lg font-semibold")
    DataTable(
        columns=[
            _col("Variant", **_VARIANT_LINK),
            _gene_column(),
            _col(
                "VEP_Consequence",
                "Consequence",
                cellType="badge_list",
                badgesField="ConsequenceBadges",
            ),
            _col(
                "VEP_CLIN_SIG",
                "ClinVar",
                cellType="badge_list",
                badgesField="ClinVarBadges",
            ),
            _col("gnomAD", "gnomAD faf95", cellType="number"),
            _col("Carriers", cellType="number"),
            _col("Families", cellType="number"),
            _col("Cohorts"),
            _col("Validations"),
            _col("Diagnostics"),
        ],
        rows=variants,
        row_key="Variant",
        pagination={"rowsPerPage": 25},
    )
    ui.label("Carriers").classes("text-lg font-semibold mt-4")
    DataTable(
        columns=[
            *_person_columns(multi),
            _col("Variant", **_VARIANT_LINK),
            _col(
                "VEP_Consequence",
                "Consequence",
                cellType="badge_list",
                badgesField="ConsequenceBadges",
            ),
            _col("GT"),
            _col("Father GT"),
            _col("Mother GT"),
            _col("Configs", "Wombat configs"),
            *_curation_columns(),
        ],
        rows=carriers,
        row_key="_key",
        pagination={"rowsPerPage": 50},
    )


def _sv_table(rows: List[Dict[str, Any]], multi: bool) -> None:
    if not rows:
        _empty("No WisecondorX call over this gene in any cohort.")
        return
    DataTable(
        columns=[
            *_person_columns(multi),
            _col("Locus"),
            _col("Type"),
            _col("wisecondorX", "Call", cellType="cnv_call"),
            _col("Size (kb)", cellType="number"),
            _col("Ratio", cellType="number"),
            _col("Exonic"),
            *_curation_columns(),
        ],
        rows=rows,
        row_key="_key",
        pagination={"rowsPerPage": 50},
    )


def _str_table(rows: List[Dict[str, Any]], genotyped: int, multi: bool) -> None:
    if not genotyped:
        _empty("No ExpansionHunter locus of this gene in any cohort catalog.")
        return
    ui.label(
        f"{len(rows)} genotype{'s' * (len(rows) != 1)} past normal, of "
        f"{genotyped} individual{'s' * (genotyped != 1)} genotyped."
    ).classes("text-sm text-gray-600")
    if not rows:
        return
    DataTable(
        columns=[
            *_person_columns(multi),
            _col("Locus"),
            _col("Alleles"),
            _col("status", "Status", cellType="badge_list", badgesField="StatusBadges"),
            _col(
                "severity",
                "Severity",
                cellType="badge_list",
                badgesField="SeverityBadges",
            ),
            _col("Disease"),
            _col("Catalog"),
            *_curation_columns(),
        ],
        rows=rows,
        row_key="_key",
        pagination={"rowsPerPage": 50},
    )


def _roh_table(rows: List[Dict[str, Any]], multi: bool) -> None:
    if not rows:
        _empty("No run of homozygosity over this gene in any cohort.")
        return
    DataTable(
        columns=[
            *_person_columns(multi),
            _col("Segment"),
            _col("Length (Mb)", cellType="number"),
            _col("Class"),
            _col("Gene covered", "Gene covered (%)", cellType="number"),
        ],
        rows=rows,
        row_key="_key",
        pagination={"rowsPerPage": 50},
    )


def _diagnostic_table(rows: List[Dict[str, Any]], multi: bool) -> None:
    if not rows:
        _empty("No diagnostic names this gene.")
        return
    DataTable(
        columns=[
            *([_col("Data dir")] if multi else []),
            _col("Type"),
            _col("FID", "Family", **_FAMILY_LINK),
            _col("Sample"),
            _col("Variant"),
            _col("Impact", cellType="badge_list", badgesField="ImpactBadges"),
            _col("Diagnostic", cellType="badge_list", badgesField="DiagnosticBadges"),
            _col("Users"),
            _col("Timestamp"),
        ],
        rows=rows,
        row_key="_key",
        pagination={"rowsPerPage": 50},
    )


# ---------------------------------------------------------------------------
# Variant page
# ---------------------------------------------------------------------------


def _load_variant(query: gene_index.VariantQuery, scorer: Any) -> Dict[str, Any]:
    dirs = data_dir_infos()
    curation = load_curation(dirs)
    carriers = snv_carrier_rows(
        gene_index.lookup_variant(query), curation, dirs, scorer
    )
    # Validations and diagnostics are keyed by the table's own spelling of the
    # chromosome: try the typed one too, for data written without "chr".
    spellings = {query.label, query.label.removeprefix("chr")}
    spellings |= {r["Variant"] for r in carriers}
    validations: List[Dict[str, Any]] = []
    diagnostics: List[Dict[str, Any]] = []
    for key in sorted(spellings):
        validations += validation_history(dirs, dirs, "SNV", key)
        diagnostics += diagnostic_rows(dirs, dirs, scorer, variant=key)
    return {
        "multi_dir": len(dirs) > 1,
        "carriers": carriers,
        "validations": validations,
        "diagnostics": diagnostics,
    }


@ui.page("/variant/{query}")
def variant_page(query: str) -> Optional[RedirectResponse]:
    """Have we seen this variant? Its carriers and curation, everywhere."""
    if redirect := check_auth():
        return redirect
    create_header()
    parsed = gene_index.parse_variant(query)
    scorer = get_gene_scorer()

    with ui.column().classes("w-full px-6 py-6 gap-3"):
        if parsed is None:
            ui.label(f"Not a variant: {query}").classes("text-xl text-red-600")
            ui.label(
                "Write it as chromosome, position, reference and alternate alleles, "
                "e.g. chr22:50674000:G:A or 22-50674000-G-A."
            ).classes("text-gray-600")
            return None
        ui.label(f"🔎 {parsed.label}").classes("text-3xl font-bold text-blue-900")
        content = ui.column().classes("w-full gap-3")

        async def load() -> None:
            content.clear()
            with content:
                ui.spinner(size="lg")
            data = await asyncio.to_thread(_load_variant, parsed, scorer)
            content.clear()
            with content:
                _render_variant(data)

        _index_line(lambda: ui.timer(0, load, once=True))
        ui.timer(0.1, load, once=True)
    return None


def _render_variant(data: Dict[str, Any]) -> None:
    carriers = data["carriers"]
    multi = data["multi_dir"]
    if carriers:
        genes = sorted({g for r in carriers for g in r["Genes"].split(", ") if g})
        with ui.row().classes("items-center gap-2"):
            ui.label("Seen in " + _seen_text(seen_summary(carriers)) + ".").classes(
                "text-base text-gray-800"
            ).mark("variant-summary")
            for gene in genes:
                ui.link(gene, f"/gene/{gene}").classes("text-sm")
        DataTable(
            columns=[
                *_person_columns(multi),
                _gene_column(),
                _col(
                    "VEP_Consequence",
                    "Consequence",
                    cellType="badge_list",
                    badgesField="ConsequenceBadges",
                ),
                _col("GT"),
                _col("Father GT"),
                _col("Mother GT"),
                _col("gnomAD", "gnomAD faf95", cellType="number"),
                _col("Configs", "Wombat configs"),
                *_curation_columns(),
            ],
            rows=carriers,
            row_key="_key",
            pagination={"rowsPerPage": 50},
        )
    else:
        ui.label("Not seen in any cohort Wombat table.").classes(
            "text-base text-gray-800"
        ).mark("variant-summary")

    ui.label("Validations").classes("text-lg font-semibold mt-4")
    if data["validations"]:
        DataTable(
            columns=[
                *([_col("Data dir")] if multi else []),
                _col("FID", "Family", **_FAMILY_LINK),
                _col("Sample"),
                _col("Validation"),
                _col("Inheritance"),
                _col("User"),
                _col("Comment"),
                _col("Ignored"),
                _col("Timestamp"),
            ],
            rows=data["validations"],
            row_key="_key",
            pagination={"rowsPerPage": 25},
        )
    else:
        _empty("No validation recorded for this variant.")

    ui.label("Diagnostics").classes("text-lg font-semibold mt-4")
    if data["diagnostics"]:
        _diagnostic_table(data["diagnostics"], multi)
    else:
        _empty("No diagnostic recorded for this variant.")
