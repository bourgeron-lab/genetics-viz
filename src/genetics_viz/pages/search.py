"""Search page for cohort-wide variant search.

The page is three parts: the form (``search_form.py``) edits the source
components and individual filters; :func:`_run_search` turns them into rows
with the logic in ``utils/search_rows.py``; the results (``search_results.py``)
show the table.
"""

import asyncio
from typing import Any, Dict, List

import polars as pl
from fastapi import Request
from nicegui import context, ui

from genetics_viz.components.filters import notify_unknown_terms
from genetics_viz.components.header import create_header
from genetics_viz.components.page_link import PageLink
from genetics_viz.pages.search_form import SearchForm, render_search_form
from genetics_viz.pages.search_results import SearchResults, render_results
from genetics_viz.utils.auth import check_auth
from genetics_viz.utils.column_names import reorder_columns_by_group
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.genesets import geneset_sources
from genetics_viz.utils.genetrek import get_snapshot
from genetics_viz.utils.inheritance import annotate_moi
from genetics_viz.utils.locus import (
    filter_bed_dataframe,
    filter_dataframe,
    parse_locus_query,
)
from genetics_viz.utils.sample_counts import MissingHeaderError
from genetics_viz.utils.search_rows import (
    DEFAULT_VISIBLE,
    PED_MISSING,
    CurationMaps,
    apply_geneset_filter,
    build_wcx_rows,
    build_wombat_rows,
    combine_wcx,
    combine_wombat,
    filter_individuals,
    filter_wcx_calls,
    find_source,
    find_sources,
    form_from_link,
    form_to_link,
    initial_preset,
    load_wombat_source,
    missing_query,
    moi_pedigree,
    pedigree_data_from_cohort,
    result_columns,
)
from genetics_viz.utils.view_presets import VIEW_PRESETS
from genetics_viz.utils.wisecondorx import parse_wisecondorx_bed


@ui.page("/search/{cohort_name}")
def search_cohort_page(cohort_name: str, request: Request) -> None:
    """Search page for cohort-wide variant search.

    Running a search writes the form into the URL query (``q``), and a page
    opened with one fills the form and runs it: the address is a link to the
    results.
    """
    if redirect := check_auth():
        return redirect
    link = PageLink.from_request(request)
    create_header(cohort_name)

    # Add IGV.js library
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.13/dist/igv.min.js"></script>
    """)
    ui.add_css("""
        .exclude-input-wrap .q-field__control:before,
        .exclude-input-wrap .q-field__control:after { border: none !important; }
    """)

    try:
        store = get_data_store()

        with ui.column().classes("w-full px-6 py-6"):
            # Title
            ui.label(f"🔍 Search: {cohort_name}").classes(
                "text-3xl font-bold text-blue-900 mb-6"
            )

            # Scan for source files (wombat + SVS)
            source_files = find_sources(store.data_dir, cohort_name)
            if not source_files:
                ui.label(
                    "No source files found (checked wombat/ and svs/ directories)"
                ).classes("text-gray-500 text-lg italic")
                return

            # Load pedigree data from the already-parsed Cohort object
            pedigree_data = pedigree_data_from_cohort(store.get_cohort(cohort_name))

            form = SearchForm(
                wombat_sources=[
                    sf for sf in source_files if sf["source_type"] == "wombat"
                ],
                wcx_sources=[
                    sf for sf in source_files if sf["source_type"] == "wisecondorx"
                ],
                # Gene sets for the filter menus: the user's tiers, GeneTrek's
                # lists and the data directory's params/genesets/*.tsv
                genesets=geneset_sources(store.data_dir, get_gene_scorer()),
                # Unique sex and phenotype values for the individual filters,
                # without the missing-value sentinels ("", "0", "-9")
                sex_values=sorted(
                    {v.get("Sex", "") for v in pedigree_data.values()} - PED_MISSING
                ),
                phenotype_values=sorted(
                    {v.get("Phenotype", "") for v in pedigree_data.values()}
                    - PED_MISSING
                ),
            )
            # The linked search, else one default wombat component
            linked = form_from_link(
                link.get("q"), form.wombat_sources, form.wcx_sources
            )
            if linked is not None and linked[0]:
                form.components[:], form.individuals = linked
            else:
                linked = None
                form.reset_components()

            search_button = render_search_form(form)

            # Results container
            results_container = ui.column().classes("w-full")

            # Capture client context for callbacks
            page_client = context.client

            async def perform_search():
                link.set(
                    "q", form_to_link(form.components, form.individuals, source_files)
                )
                link.push()
                await _run_search(
                    form,
                    results_container,
                    store=store,
                    cohort_name=cohort_name,
                    source_files=source_files,
                    pedigree_data=pedigree_data,
                    page_client=page_client,
                )

            # Set up handlers after function definition
            search_button.on_click(perform_search)
            form.set_search(perform_search)
            if linked is not None:
                ui.timer(0.1, perform_search, once=True)

    except Exception as e:
        with ui.column().classes("w-full px-6 py-6"):
            _show_error(e)


def _show_error(e: Exception) -> None:
    import traceback

    ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
    ui.label("Traceback:").classes("text-red-500 font-semibold")
    ui.label(traceback.format_exc()).classes(
        "text-red-500 text-xs font-mono whitespace-pre"
    )


def _show_message(container: Any, text: str, classes: str) -> None:
    container.clear()
    with container:
        ui.label(text).classes(classes)


async def _run_search(
    form: SearchForm,
    results_container: Any,
    *,
    store: Any,
    cohort_name: str,
    source_files: List[Dict[str, Any]],
    pedigree_data: Dict[str, Dict[str, str]],
    page_client: Any,
) -> None:
    """Execute the search — iterate over source components, union results."""
    results_container.clear()
    components = form.components

    # Show progress indicator while loading
    with results_container:
        with ui.column().classes("items-center gap-4 justify-center py-8 w-full"):
            progress = ui.circular_progress(
                min=0, max=100, value=0, size="xl", color="blue"
            )
            status_label = ui.label("Starting search...").classes(
                "text-lg text-gray-600"
            )

    # Validate: at least one component, each with a locus or genesets
    problem = missing_query(components)
    if problem:
        _show_message(results_container, problem, "text-orange-600")
        return

    async def step(value: int, text: str) -> None:
        progress.set_value(value)
        status_label.set_text(text)
        await asyncio.sleep(0)

    await step(10, "Loading data...")

    try:
        # --- Process each component ---
        all_wombat_dfs: List[pl.DataFrame] = []
        all_wcx_dfs: List[pl.DataFrame] = []
        cohort = store.get_cohort(cohort_name)

        for comp in components:
            sf = find_source(source_files, comp["source_key"])
            if sf is None:
                continue

            locus_val = comp["locus"]["value"]
            query_params = None
            if locus_val and locus_val.strip():
                query_params = parse_locus_query(locus_val)

            if comp["type"] == "wombat":
                try:
                    df = await asyncio.to_thread(
                        load_wombat_source, store.data_dir, cohort, sf
                    )
                except MissingHeaderError as exc:
                    # Read anyway, its first row would be the columns
                    results_container.clear()
                    with results_container:
                        ui.label(
                            f"Cannot search Wombat source {sf['display_name']}: {exc}"
                        ).classes("text-red-600").mark("search-source-error")
                    return
                if df is None or len(df) == 0:
                    continue
                # Locus filter
                if query_params:
                    df = filter_dataframe(df, query_params)
                # Geneset filter
                df = apply_geneset_filter(df, comp, "wombat", form.genesets)
                if df is not None and len(df) > 0:
                    all_wombat_dfs.append(df)

            elif comp["type"] == "wisecondorx":
                df = await asyncio.to_thread(parse_wisecondorx_bed, sf["file_path"])
                if df is None or len(df) == 0:
                    continue
                # Locus filter (use exonic columns if exonic_only is checked)
                exonic = comp.get("exonic_only", {}).get("value", False)
                if query_params:
                    df = filter_bed_dataframe(df, query_params, exonic=exonic)
                # Geneset filter, then calls and |ratio|
                df = apply_geneset_filter(df, comp, "wisecondorx", form.genesets)
                df = filter_wcx_calls(df, comp)
                if df is not None and len(df) > 0:
                    all_wcx_dfs.append(df)

        await step(35, "Filtering data...")

        # Deduplicate rows across components
        wombat_filtered = combine_wombat(all_wombat_dfs)
        wcx_filtered = combine_wcx(all_wcx_dfs)

        await step(60, "Processing variants...")

        wombat_count = len(wombat_filtered) if wombat_filtered is not None else 0
        wcx_count = len(wcx_filtered) if wcx_filtered is not None else 0

        if wombat_count == 0 and wcx_count == 0:
            _show_message(
                results_container, "No results found", "text-gray-500 text-lg italic"
            )
            return

        # Load validation and diagnostic data
        maps = CurationMaps(store.data_dir).load(
            snv=wombat_filtered is not None, sv=wcx_filtered is not None
        )

        # Yield to event loop before badge processing
        await asyncio.sleep(0)

        # Track unknown terms for warnings
        unknown_consequences: set = set()
        unknown_clinvar_terms: set = set()
        sample_to_family = {sid: ped["FID"] for sid, ped in pedigree_data.items()}
        common: Dict[str, Any] = {
            "cohort_name": cohort_name,
            "pedigree_data": pedigree_data,
            "sample_to_family": sample_to_family,
            "maps": maps,
            "gene_scorer": get_gene_scorer(),
        }

        # Each source's rows are filtered by its own components: a row passes
        # if ANY component of its type would keep it.
        all_rows: List[Dict[str, Any]] = []
        if wombat_filtered is not None and wombat_count > 0:
            all_rows.extend(
                build_wombat_rows(
                    wombat_filtered,
                    [c for c in components if c["type"] == "wombat"],
                    unknown_consequences=unknown_consequences,
                    unknown_clinvar=unknown_clinvar_terms,
                    **common,
                )
            )
        if wcx_filtered is not None and wcx_count > 0:
            all_rows.extend(
                build_wcx_rows(
                    wcx_filtered,
                    [c for c in components if c["type"] == "wisecondorx"],
                    **common,
                )
            )

        await step(85, "Rendering table...")

        notify_unknown_terms(unknown_consequences, unknown_clinvar_terms)

        # Group same-group columns together
        all_columns = reorder_columns_by_group(
            result_columns(wombat_filtered, wcx_filtered)
        )

        # Auto-select preset based on component mix; its columns, else defaults
        preset = initial_preset(components, source_files, VIEW_PRESETS)
        initial_selected = [
            col for col in preset.get("columns", []) if col in all_columns
        ]
        selected_cols = {
            "value": initial_selected
            or [col for col in DEFAULT_VISIBLE if col in all_columns]
        }

        all_rows = filter_individuals(all_rows, pedigree_data, form.individuals)

        # MOI: does each SNV fit its gene's mode of inheritance?
        # Second hits are looked for among these results only.
        annotate_moi(all_rows, get_snapshot(), pedigree=moi_pedigree(pedigree_data))

        await step(100, "Complete!")

        # Clear progress indicator and show results
        results_container.clear()
        with results_container:
            render_results(
                SearchResults(
                    rows=all_rows,
                    all_columns=all_columns,
                    selected_cols=selected_cols,
                    selected_preset={"name": preset["name"]},
                    has_wcx=any(c["type"] == "wisecondorx" for c in components),
                    maps=maps,
                    cohort_name=cohort_name,
                    sample_to_family=sample_to_family,
                    data_dir=store.data_dir,
                    page_client=page_client,
                )
            )

    except Exception as e:
        results_container.clear()
        with results_container:
            _show_error(e)
