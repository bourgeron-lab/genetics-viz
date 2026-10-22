"""Shared component for SV visualization in a dialog with IGV.js.

``show_sv_dialog`` builds the dialog and wires its panels together:

- the samples shown (carrier, parents, others added);
- SV details and suggested positions (``sv_details.py``);
- Read Evidence (``sv_read_panel.py``), then the Copy Number View
  (``sv_cn_panel.py``) or its bedgraph fallback, the VAF and read-level views
  (``sv_igv_views.py``);
- the validation form (``sv_validation_panel.py``).

The views share one ``SVViewState``. Its ``bounds`` dict is the SV's current
position: every curated change moves it, and the read views follow.
"""

from typing import Any, Dict, List, Optional, Tuple

from nicegui import ui

from genetics_viz.components.sv_cn_panel import CopyNumberPanel
from genetics_viz.components.sv_details import (
    build_sv_suggestions,
    render_suggestions,
    render_sv_details,
)
from genetics_viz.components.sv_igv_views import (
    ReadLevelView,
    SVViewState,
    coverage_view,
    has_bedgraph,
    update_roi_js,
    vaf_view,
)
from genetics_viz.components.sv_read_panel import ReadEvidencePanel
from genetics_viz.components.sv_validation_panel import (
    SVValidationForm,
    render_sv_validation,
)
from genetics_viz.utils.cnv_coverage import CNRequest, find_sample_coverage
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.pipeline_params import load_params
from genetics_viz.utils.wisecondorx import infer_sv_type

_PARENT_MISSING = ("-", "0")


def _expanded_locus(chrom: str, start: int, end: int) -> str:
    """The locus 1.6 times as long as the SV, centred on it."""
    center = (start + end) // 2
    expanded_length = int((end - start) * 1.6)
    expanded_start = max(0, center - (expanded_length // 2))
    return f"{chrom}:{expanded_start}-{center + (expanded_length // 2)}"


def _curated_position(
    validations: List[Tuple],
) -> Tuple[Optional[str], Optional[str]]:
    """Curated start and end of the newest present validation that has one.

    Keep original coordinates for the title and validation lookups; curated
    ones are only for the ROI and the IGV locus.
    """
    present_with_curated = [
        v for v in validations if v[0] == "present" and v[3] != "1" and (v[4] or v[5])
    ]
    if not present_with_curated:
        return None, None
    present_with_curated.sort(key=lambda v: v[6], reverse=True)
    most_recent = present_with_curated[0]
    return most_recent[4] or None, most_recent[5] or None


def show_sv_dialog(
    cohort_name: str,
    family_id: str,
    chrom: str,
    start: str,
    end: str,
    sample: str,
    sv_data: Dict[str, Any],
    on_validation_saved: Optional[callable] = None,
    family_members_override: Optional[List[str]] = None,
    sample_parents_override: Optional[Dict[str, Optional[str]]] = None,
) -> None:
    """Show SV visualization dialog with IGV viewer.

    Args:
        cohort_name: Cohort name
        family_id: Family ID
        chrom: Chromosome
        start: Start position
        end: End position
        sample: Sample ID
        sv_data: Additional SV data to display
        on_validation_saved: Optional callback to invoke after validation is saved
        family_members_override: Optional pre-computed family member IDs (for standalone pages)
        sample_parents_override: Optional pre-computed parent IDs (for standalone pages)
    """
    store = get_data_store()
    start_pos = int(start)
    end_pos = int(end)
    locus = f"{chrom}:{start}-{end}"

    # Check for curated boundaries from validations
    from genetics_viz.components.validation_loader import load_validation_map

    validation_file = store.data_dir / "validations" / "svs.tsv"
    validation_map = load_validation_map(validation_file, family_id)
    all_validation_map = load_validation_map(validation_file, family_id=None)

    # Determine SV type for variant key (dup for gain, del for loss)
    sv_type = infer_sv_type(sv_data)
    variant_key = f"{chrom}:{start}-{end}:{sv_type}"

    curated_start, curated_end = _curated_position(
        validation_map.get((variant_key, sample), [])
    )
    curated_locus = None
    igv_locus = _expanded_locus(chrom, start_pos, end_pos)
    if curated_start or curated_end:
        display_start = curated_start or start
        display_end = curated_end or end
        curated_locus = f"{chrom}:{display_start}-{display_end}"
        igv_locus = _expanded_locus(chrom, int(display_start), int(display_end))

    with (
        ui.dialog().props("maximized") as dialog,
        ui.card().classes("w-full h-full"),
    ):
        with ui.column().classes("w-full h-full p-6"):
            # Header with close button
            with ui.row().classes("items-center justify-between w-full mb-4"):
                with ui.row().classes("items-center gap-2"):
                    ui.label(f"🧬 SV: {locus} - {sample}").classes(
                        "text-2xl font-bold text-blue-900"
                    )
                    if curated_locus:
                        ui.label(f"(curated: {curated_locus})").classes(
                            "text-lg font-semibold text-green-600"
                        )
                ui.button(icon="close", on_click=lambda: dialog.close()).props(
                    "flat round"
                )

            # Get family members (use overrides if provided, else query cohort)
            if family_members_override is not None:
                family_members = family_members_override
            else:
                family_members = []
            if sample_parents_override is not None:
                sample_parents = sample_parents_override
            else:
                sample_parents = {"father": None, "mother": None}

            cohort = store.get_cohort(cohort_name)
            if family_members_override is None or sample_parents_override is None:
                if cohort:
                    members_data = cohort.get_family_members(family_id)
                    if family_members_override is None:
                        family_members = [m["Sample ID"] for m in members_data]
                    if sample_parents_override is None:
                        for member in members_data:
                            if member["Sample ID"] == sample:
                                sample_parents["father"] = member.get("Father")
                                sample_parents["mother"] = member.get("Mother")
                                break

            def has_coverage(sample_id: str) -> bool:
                """A coverage track can be drawn: bigWig or by1000 bedgraph."""
                return has_bedgraph(store.data_dir, sample_id) or bool(
                    find_sample_coverage(store.data_dir, sample_id)
                )

            def cn_request(sample_id: str) -> CNRequest:
                """The sample's family and pedigree sex, for its panel group."""
                if cohort:
                    for fam in cohort.families.values():
                        member = fam.get_sample(sample_id)
                        if member is not None:
                            return CNRequest(sample_id, fam.family_id, member.sex)
                return CNRequest(sample_id, family_id, None)

            def get_relationship_label(sample_id: str) -> str:
                """Get relationship label for a sample."""
                if sample_id == sample:
                    return "(carrier)"

                if sample_id == sample_parents["father"]:
                    return "(father)"
                if sample_id == sample_parents["mother"]:
                    return "(mother)"

                # Check if it's a sibling
                if cohort and sample_parents["father"] and sample_parents["mother"]:
                    members_data = cohort.get_family_members(family_id)
                    for member in members_data:
                        if member["Sample ID"] == sample_id:
                            member_father = member.get("Father")
                            member_mother = member.get("Mother")
                            if (
                                member_father == sample_parents["father"]
                                and member_mother == sample_parents["mother"]
                                and member_father
                                and member_mother
                                and member_father not in _PARENT_MISSING
                                and member_mother not in _PARENT_MISSING
                            ):
                                return "(sibling)"
                            break

                return ""

            state = SVViewState(
                data_dir=store.data_dir,
                chrom=chrom,
                start=start,
                end=end,
                sample=sample,
                sv_type=sv_type,
                call=sv_data.get("call", ""),
                igv_locus=igv_locus,
                roi_start=curated_start or start,
                roi_end=curated_end or end,
                label_for=get_relationship_label,
            )

            # Automatically add parents if available
            for parent_id in sample_parents.values():
                if (
                    parent_id
                    and parent_id not in _PARENT_MISSING
                    and parent_id != sample
                    and has_coverage(parent_id)
                ):
                    state.additional.append(parent_id)

            # Additional samples section with add menu
            with ui.row().classes("items-center gap-4 mb-2"):
                # Menu to add samples
                with ui.button("Add Samples", icon="add").props(
                    "outline color=blue size=sm"
                ):
                    with ui.menu():
                        ui.menu_item("Add Parents", on_click=lambda: add_parents())
                        ui.menu_item("Add Family", on_click=lambda: add_family())
                        ui.separator()
                        with ui.row().classes("items-center gap-2 px-4 py-2"):
                            barcode_input = (
                                ui.input("Barcode").classes("flex-grow").props("dense")
                            )
                            ui.button(
                                "Add",
                                icon="add",
                                on_click=lambda: add_sample(barcode_input.value),
                            ).props("flat dense size=sm")

            # Display additional samples
            additional_samples_container = ui.column().classes("gap-1 mb-4")

            def refresh_additional_samples():
                additional_samples_container.clear()
                with additional_samples_container:
                    if state.additional:
                        ui.label("Additional Samples:").classes(
                            "text-sm font-semibold text-gray-700"
                        )
                        for add_sample_id in state.additional:
                            with ui.row().classes("items-center gap-2"):
                                ui.label(state.label(add_sample_id)).classes(
                                    "text-sm text-gray-600"
                                )
                                ui.button(
                                    icon="delete",
                                    on_click=lambda _, sid=add_sample_id: remove_sample(
                                        sid
                                    ),
                                ).props("flat dense size=xs color=red")

            def samples_changed():
                refresh_additional_samples()
                refresh_igv()
                vaf.refresh()
                read_level.refresh()
                refresh_read_evidence()

            def add_new(candidates: List[str]) -> List[str]:
                """Add the candidates not shown yet that have coverage."""
                added = []
                for sid in candidates:
                    if sid not in state.additional and sid != sample:
                        if has_coverage(sid):
                            state.additional.append(sid)
                            added.append(sid)
                return added

            def add_sample(sample_id: str):
                if (
                    sample_id
                    and sample_id not in state.additional
                    and sample_id != sample
                ):
                    if add_new([sample_id]):
                        samples_changed()
                    else:
                        ui.notify(
                            f"No coverage file found for sample: {sample_id}",
                            type="warning",
                        )

            def add_parents():
                added = add_new(
                    [
                        p
                        for p in sample_parents.values()
                        if p and p not in _PARENT_MISSING
                    ]
                )
                if added:
                    samples_changed()
                    ui.notify(f"Added parents: {', '.join(added)}", type="positive")
                else:
                    ui.notify("No parents to add or files not found", type="warning")

            def add_family():
                added = add_new(family_members)
                if added:
                    samples_changed()
                    ui.notify(f"Added {len(added)} family members", type="positive")
                else:
                    ui.notify("No additional family members to add", type="warning")

            def remove_sample(sample_id: str):
                if sample_id in state.additional:
                    state.additional.remove(sample_id)
                    samples_changed()

            refresh_additional_samples()

            render_sv_details(sv_data, locus, get_gene_scorer())

            # Filled once the validation form is built, further down.
            form: Dict[str, Optional[SVValidationForm]] = {"value": None}

            def set_curated(key: str, value: int) -> None:
                state.bounds[key] = int(value)
                if form["value"] is not None:
                    form["value"].curated_input(key).value = str(value)

            def apply_suggestion(s: Dict[str, Any]) -> None:
                set_curated("start", s["start"])
                set_curated("end", s["end"])
                update_roi()
                # Preselect inheritance from parent relationship
                inh = s.get("inheritance", "unknown")
                if inh != "unknown" and form["value"] is not None:
                    form["value"].inheritance.value = inh
                ui.notify(
                    f"Curated position set to {chrom}:{s['start']}-{s['end']}",
                    type="positive",
                )

            # Suggested coordinates from overlapping validated SVs
            render_suggestions(
                build_sv_suggestions(
                    all_validation_map=all_validation_map,
                    chrom=chrom,
                    start_pos=start_pos,
                    end_pos=end_pos,
                    current_sample=sample,
                    sample_parents=sample_parents,
                    family_members=family_members,
                ),
                chrom,
                apply_suggestion,
            )

            cohort_params = (
                load_params(cohort.params_file) if cohort and cohort.params_file else {}
            )
            # Read Evidence comes first. It is built further down, once the read
            # windows are defined, into this slot.
            read_slot = ui.column().classes("w-full gap-0")

            # The carrier's normalised coverage bigWig, divided by the cohort's
            # control panel, replaces the by1000 bedgraph in IGV when present.
            # Folded then: Read Evidence already draws the copy number.
            use_cn = bool(find_sample_coverage(store.data_dir, sample))
            cn_state = {"loaded": False}
            cnv_expansion = (
                ui.expansion(
                    "Copy Number View" if use_cn else "CNV Coverage View",
                    icon="analytics",
                    value=not use_cn,
                )
                .classes("w-full mt-4 border border-gray-400 rounded-lg")
                .props("header-class='text-lg font-bold bg-gray-100'")
                .mark("cn-expansion")
            )

            with cnv_expansion:
                if use_cn:
                    cn_panel = CopyNumberPanel(
                        data_dir=store.data_dir,
                        cohort_dir=cohort.path if cohort else None,
                        cohort_name=cohort_name,
                        params=cohort_params,
                        chrom=chrom,
                        roi=(int(state.roi_start), int(state.roi_end)),
                        call=str(sv_data.get("call", "")),
                        label_for=get_relationship_label,
                    )

                    def refresh_igv():
                        """Redraw the copy-number panel for the current samples.

                        Folded, it is drawn on first opening instead: a chart
                        built in a hidden panel has no width.
                        """
                        cn_state["loaded"] = bool(cnv_expansion.value)
                        if not cn_state["loaded"]:
                            return
                        cn_panel.refresh([cn_request(sid) for sid in state.samples()])
                else:
                    coverage = coverage_view(state)
                    refresh_igv = coverage.refresh

                refresh_igv()

            # Handle expansion panel opening to resize IGV
            def on_cnv_expansion_change(e):
                if e.value and use_cn and not cn_state["loaded"]:
                    refresh_igv()
                elif e.value:  # Panel was opened
                    ui.run_javascript("""
                        setTimeout(() => {
                            if (window.igvBrowser && window.igvBrowser.visibilityChange) {
                                window.igvBrowser.visibilityChange();
                            }
                        }, 100);
                    """)

            cnv_expansion.on_value_change(on_cnv_expansion_change)

            vaf_expansion = (
                ui.expansion("VAF View", icon="water_drop", value=False)
                .classes("w-full mt-4 border border-gray-400 rounded-lg")
                .props("header-class='text-lg font-bold bg-gray-100'")
            )
            with vaf_expansion:
                vaf = vaf_view(state)
                vaf.refresh()
            vaf_expansion.on_value_change(lambda e: vaf.resize() if e.value else None)

            # Read-level view: the CRAMs in the read windows (split view)
            cram_expansion = (
                ui.expansion("Read-Level View", icon="dna")
                .classes("w-full mt-6 border border-gray-400 rounded-lg")
                .props("default-opened header-class='text-lg font-bold bg-gray-100'")
            )

            def on_pick(key: str, position: int) -> None:
                set_curated(key, position)
                update_roi()

            with cram_expansion:
                read_level = ReadLevelView(state, on_pick)
            cram_expansion.on_value_change(
                lambda e: read_level.resize() if e.value else None
            )

            def update_roi():
                """Move every view to the current bounds."""
                ui.run_javascript(update_roi_js(state))
                refresh_read_evidence()
                if use_cn:
                    cn_panel.set_roi((state.bounds["start"], state.bounds["end"]))
                    if cn_state["loaded"]:
                        refresh_igv()

            def on_read_mode_change(e):
                state.read_mode = e.value
                read_level.refresh()
                refresh_read_evidence()

            # Samplot-like view of the same windows: discordant pairs at their
            # insert size, split reads, soft clips and the copy number. Shown
            # first, above the Copy Number View (read_slot).
            with read_slot:
                reads_expansion = (
                    ui.expansion("Read Evidence", icon="insights", value=True)
                    .classes("w-full mt-4 border border-gray-400 rounded-lg")
                    .props("header-class='text-lg font-bold bg-gray-100'")
                )
                with reads_expansion:
                    with ui.row().classes("items-center gap-3 px-2 pb-2"):
                        ui.toggle(
                            {
                                "auto": "Auto",
                                "single": "Single window",
                                "split": "Split",
                            },
                            value="auto",
                            on_change=on_read_mode_change,
                        ).props("dense no-caps size=sm").mark("read-mode")
                        read_evidence_label = ui.label().classes(
                            "text-sm text-gray-500"
                        )

                        def on_low_mapq_change(e):
                            state.low_mapq = bool(e.value)
                            refresh_read_evidence()

                        ui.checkbox(
                            "Show MAPQ 0 reads (faded)",
                            value=False,
                            on_change=on_low_mapq_change,
                        ).props("dense").classes("text-sm").mark("read-low-mapq")
                    read_panel = ReadEvidencePanel(
                        data_dir=store.data_dir,
                        cohort_dir=cohort.path if cohort else None,
                        cohort_name=cohort_name,
                        params=cohort_params,
                        label_for=get_relationship_label,
                        on_use_breakpoints=lambda s_, e_, approx: use_read_breakpoints(
                            s_, e_, approx
                        ),
                    )

            def refresh_read_evidence():
                used_mode, windows = state.read_windows()
                read_evidence_label.text = state.describe_read_mode(used_mode)
                read_panel.refresh(
                    [cn_request(sid) for sid in state.samples()],
                    windows,
                    [state.bounds["start"], state.bounds["end"]],
                    include_low_mapq=state.low_mapq,
                    sv_type=sv_type,
                )

            def use_read_breakpoints(new_start, new_end, approximate=False):
                """Take the suggested breakpoints as the curated position."""
                for key, value in (("start", new_start), ("end", new_end)):
                    if value is not None:
                        set_curated(key, value)
                update_roi()
                ui.notify(
                    f"Curated position set to {chrom}:{state.bounds['start']}"
                    f"-{state.bounds['end']}"
                    + (
                        " (approximate, from the coverage: check it)"
                        if approximate
                        else " from the reads"
                    ),
                    type="warning" if approximate else "positive",
                )

            refresh_read_evidence()

            def on_saved():
                # Call the callback if provided, then close dialog
                if on_validation_saved:
                    on_validation_saved()
                dialog.close()

            form["value"] = render_sv_validation(
                validation_file=validation_file,
                family_id=family_id,
                variant_key=variant_key,
                chrom=chrom,
                start=start,
                end=end,
                sample=sample,
                sample_parents=sample_parents,
                bounds=state.bounds,
                on_bounds_changed=update_roi,
                on_saved=on_saved,
            )

        dialog.open()
