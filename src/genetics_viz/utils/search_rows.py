"""Search page logic without UI: sources, per-source filters, result rows.

The Search page (``pages/search.py``) unions the rows of several *source
components*, each a Wombat cohort table or the cohort's WisecondorX calls with
its own filters. Everything here is plain data in, plain data out, so the page
keeps only the widgets and the order of the steps.

A component is a dict of mutable ``{"value": ...}`` boxes (the filter widgets
write into them), built by :func:`new_wombat_component` /
:func:`new_wcx_component`.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple
from uuid import uuid4

import polars as pl

from genetics_viz.components.diagnostic_loader import (
    add_diagnostic_status_to_row,
    load_diagnostic_map,
)
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.utils.column_names import (
    genomic_sort_key,
    get_dropped_columns,
    get_schema_overrides,
)
from genetics_viz.utils.export import is_export_column
from genetics_viz.utils.gene_scoring import badge_labels, build_gene_badges, split_genes
from genetics_viz.utils.genesets import geneset_filter_expr, selected_symbols
from genetics_viz.utils.sample_counts import (
    cohort_wcx_file,
    cohort_wombat_sources,
    family_header_sources,
    read_cohort_wombat,
)
from genetics_viz.utils.validation_status import FILTER_STATUSES, status_selected
from genetics_viz.utils.variant_rows import (
    add_score_badges,
    add_vep_badges,
    consequences_selected,
    gnomad_filtered,
    in_gnomad_wgs,
    in_lcr,
)
from genetics_viz.utils.vep import VEP_CONSEQUENCES
from genetics_viz.utils.view_presets import select_preset_for_config
from genetics_viz.utils.wisecondorx import infer_sv_type
from genetics_viz.utils.wombat import aggregate_by_variant

# Values that represent "unknown/missing" in pedigree files
PED_MISSING = {"", "0", "-9"}

# WisecondorX call values (from svs_tab.py pattern)
ALL_CALL_VALUES = [
    "Robust LOSS",
    "Robust GAIN",
    "Permissive LOSS",
    "Permissive Gain",
    "Below threshold",
]

# Default visible columns when no preset names any
DEFAULT_VISIBLE = [
    "Variant",
    "VEP_Consequence",
    "gene",
    "MOI",
    "VEP_CLIN_SIG",
    "fafmax_faf95_max_genomes",
    "FID",
    "Phenotype",
    "sample",
    "sample_gt",
    "father_gt",
    "mother_gt",
    "Validation",
]


# ---------------------------------------------------------------------------
# Sources and components
# ---------------------------------------------------------------------------


def find_sources(data_dir: Path, cohort_name: str) -> List[Dict[str, Any]]:
    """The cohort's searchable files: Wombat tables, then WisecondorX calls.

    Wombat results are sorted, so the first (the default) is stable.
    """
    source_files: List[Dict[str, Any]] = [
        {
            "file_path": src.path,
            "vep_config": src.vep_config,
            "wombat_config": src.wombat_config,
            "display_name": src.wombat_config,
            "source_type": "wombat",
            "parent_dir": "wombat",
        }
        for src in cohort_wombat_sources(data_dir, cohort_name)
    ]

    # Scan for SVS subdirectories
    svs_dir = data_dir / "cohorts" / cohort_name / "svs"
    if svs_dir.exists():
        for subdir in sorted(svs_dir.iterdir()):
            if not subdir.is_dir():
                continue
            if subdir.name == "wisecondorx":
                bed_file = cohort_wcx_file(data_dir, cohort_name)
                if bed_file is not None:
                    source_files.append(
                        {
                            "file_path": bed_file,
                            "vep_config": "",
                            "wombat_config": "",
                            "display_name": subdir.name,
                            "source_type": "wisecondorx",
                            "parent_dir": "svs",
                        }
                    )
    return source_files


def source_key(sf: Dict[str, Any]) -> str:
    """The key a component's source select stores."""
    return f"{sf['parent_dir']}/{sf['display_name']}"


def find_source(
    source_files: Iterable[Dict[str, Any]], key: str
) -> Optional[Dict[str, Any]]:
    """The source file of a component's ``source_key``, or ``None``."""
    return next((s for s in source_files if source_key(s) == key), None)


def new_wombat_component(source: str = "") -> Dict[str, Any]:
    """Create a new wombat source filter component state dict."""
    return {
        "id": str(uuid4()),
        "type": "wombat",
        "source_key": source,
        "locus": {"value": ""},
        "genesets": {"value": []},
        "impacts": {"value": list(VEP_CONSEQUENCES.keys())},
        "validations": {"value": list(FILTER_STATUSES)},
        "exclude_lcr": {"value": True},
        "exclude_gnomad": {"value": True},
        "exclude_gnomad_wgs": {"value": False},
    }


def new_wcx_component(source: str = "") -> Dict[str, Any]:
    """Create a new wisecondorx source filter component state dict."""
    return {
        "id": str(uuid4()),
        "type": "wisecondorx",
        "source_key": source,
        "locus": {"value": ""},
        "genesets": {"value": []},
        "selected_calls": {
            "value": [c for c in ALL_CALL_VALUES if c != "Below threshold"]
        },
        "exonic_only": {"value": False},
        "validations": {"value": list(FILTER_STATUSES)},
        "ratio_min": {"value": None},
        "ratio_max": {"value": None},
    }


def missing_query(components: List[Dict[str, Any]]) -> Optional[str]:
    """Why the search cannot run, or ``None``.

    Every component needs a locus or a geneset: a whole cohort table would be
    too large to show.
    """
    if not components:
        return "Please add at least one source component"
    for i, comp in enumerate(components):
        locus_val = comp["locus"]["value"]
        has_locus = bool(locus_val and locus_val.strip())
        if not has_locus and not comp["genesets"]["value"]:
            return (
                f"Component {i + 1} ({comp['type']}): please enter a locus or "
                "select at least one geneset"
            )
    return None


# ---------------------------------------------------------------------------
# Loading and frame-level filters
# ---------------------------------------------------------------------------


def load_wombat_source(
    data_dir: Path, cohort: Any, sf: Dict[str, Any]
) -> Optional[pl.DataFrame]:
    """Read a cohort Wombat table (blocking; run it off-thread).

    A headerless cohort table borrows a family's header; with none to borrow
    this raises ``MissingHeaderError``.
    """
    df = read_cohort_wombat(
        sf["file_path"],
        family_header_sources(data_dir, cohort, sf["vep_config"], sf["wombat_config"])
        if cohort is not None
        else (),
        infer_schema_length=10000,
        schema_overrides=get_schema_overrides(),
        null_values=[".", ""],
    )
    if df is None:
        return None
    drop = get_dropped_columns() & set(df.columns)
    if drop:
        df = df.drop(list(drop))
    return df


def apply_geneset_filter(
    df: pl.DataFrame, comp: Dict[str, Any], source_type: str, genesets: List[Any]
) -> pl.DataFrame:
    """Keep the rows on a gene of the component's selected gene sets."""
    if not comp["genesets"]["value"]:
        return df
    wanted = selected_symbols(genesets, comp["genesets"]["value"])
    if not wanted:
        return df
    if source_type == "wombat":
        sym_col = "VEP_SYMBOL"
    else:
        exonic = comp.get("exonic_only", {}).get("value", False)
        sym_col = "exonic_symbol" if exonic else "genic_symbol"
    if sym_col not in df.columns:
        return df
    return df.filter(geneset_filter_expr(sym_col, wanted))


def filter_wcx_calls(df: pl.DataFrame, comp: Dict[str, Any]) -> pl.DataFrame:
    """The component's call types and ``|ratio|`` bounds."""
    if "wisecondorX" in df.columns and comp.get("selected_calls"):
        selected_calls = comp["selected_calls"]["value"]
        if selected_calls:
            df = df.filter(pl.col("wisecondorX").is_in(selected_calls))
    if "ratio" in df.columns:
        ratio = pl.col("ratio").cast(pl.Float64, strict=False).abs()
        ratio_min = comp.get("ratio_min", {}).get("value")
        ratio_max = comp.get("ratio_max", {}).get("value")
        if ratio_min is not None:
            df = df.filter(ratio >= ratio_min)
        if ratio_max is not None:
            df = df.filter(ratio <= ratio_max)
    return df


def combine_wombat(frames: List[pl.DataFrame]) -> Optional[pl.DataFrame]:
    """The components' Wombat rows, one per variant and sample."""
    if not frames:
        return None
    return aggregate_by_variant(pl.concat(frames, how="diagonal_relaxed"))


def combine_wcx(frames: List[pl.DataFrame]) -> Optional[pl.DataFrame]:
    """The components' WisecondorX rows, one per call and sample."""
    if not frames:
        return None
    combined = pl.concat(frames, how="diagonal_relaxed")
    dedup_cols = [
        c for c in ["chr", "start", "end", "type", "sample"] if c in combined.columns
    ]
    if dedup_cols:
        combined = combined.unique(subset=dedup_cols, keep="first", maintain_order=True)
    return combined


# ---------------------------------------------------------------------------
# Pedigree
# ---------------------------------------------------------------------------


def pedigree_data_from_cohort(cohort: Any) -> Dict[str, Dict[str, str]]:
    """Build pedigree lookup from the already-parsed Cohort object.

    Uses the DataStore's Cohort (which handles all header formats robustly)
    instead of re-parsing the pedigree file.

    Returns dict mapping sample_id -> {FID, Father, Mother, Sex, Phenotype}.
    """
    if cohort is None:
        return {}

    pedigree_data: Dict[str, Dict[str, str]] = {}
    for family in cohort.families.values():
        for sample in family.samples:
            pedigree_data[sample.sample_id] = {
                "FID": sample.family_id,
                "Father": sample.father_id or "",
                "Mother": sample.mother_id or "",
                "Sex": sample.sex or "",
                "Phenotype": sample.phenotype or "",
            }
    return pedigree_data


def moi_pedigree(pedigree_data: Dict[str, Dict[str, str]]):
    """Pedigree lookup for ``annotate_moi`` over ``pedigree_data_from_cohort``."""

    def lookup(sample_id: str) -> Dict[str, Optional[str]]:
        ped = pedigree_data.get(sample_id)
        if ped is None:
            return {}
        father = pedigree_data.get(ped.get("Father", ""), {})
        mother = pedigree_data.get(ped.get("Mother", ""), {})
        return {
            "sex": ped.get("Sex"),
            "father_pheno": father.get("Phenotype"),
            "mother_pheno": mother.get("Phenotype"),
        }

    return lookup


@dataclass
class IndividualFilters:
    """The Individuals panel: one set for every component."""

    sex: Dict[str, List[str]] = field(default_factory=lambda: {"value": []})
    phenotype: Dict[str, List[str]] = field(default_factory=lambda: {"value": []})
    has_parents: Dict[str, bool] = field(default_factory=lambda: {"value": False})
    exclude_samples: Dict[str, List[str]] = field(default_factory=lambda: {"value": []})

    def clear(self) -> None:
        self.sex["value"] = []
        self.phenotype["value"] = []
        self.has_parents["value"] = False
        self.exclude_samples["value"] = []


def filter_individuals(
    rows: List[Dict[str, Any]],
    pedigree_data: Dict[str, Dict[str, str]],
    filters: IndividualFilters,
) -> List[Dict[str, Any]]:
    """Apply individual filters (sex, phenotype, has-parents, excluded samples)."""

    def ped(r: Dict[str, Any]) -> Dict[str, str]:
        return pedigree_data.get(r.get("sample", ""), {})

    if filters.sex["value"]:
        rows = [r for r in rows if ped(r).get("Sex", "") in filters.sex["value"]]

    if filters.phenotype["value"]:
        rows = [
            r for r in rows if ped(r).get("Phenotype", "") in filters.phenotype["value"]
        ]

    if filters.has_parents["value"]:
        rows = [
            r
            for r in rows
            if ped(r).get("Father", "") not in PED_MISSING
            and ped(r).get("Mother", "") not in PED_MISSING
        ]

    if filters.exclude_samples["value"]:
        exclude_set = set(filters.exclude_samples["value"])
        rows = [r for r in rows if r.get("sample", "") not in exclude_set]
    return rows


# ---------------------------------------------------------------------------
# Curation lookups and SV keys
# ---------------------------------------------------------------------------


@dataclass
class CurationMaps:
    """Validation and diagnostic lookups of the whole data directory."""

    data_dir: Path
    snv_validation: Dict = field(default_factory=dict)
    sv_validation: Dict = field(default_factory=dict)
    snv_diagnostic: Dict = field(default_factory=dict)
    sv_diagnostic: Dict = field(default_factory=dict)

    @property
    def snv_validation_file(self) -> Path:
        return self.data_dir / "validations" / "snvs.tsv"

    @property
    def sv_validation_file(self) -> Path:
        return self.data_dir / "validations" / "svs.tsv"

    @property
    def snv_diagnostic_file(self) -> Path:
        return self.data_dir / "diagnostics" / "snvs.tsv"

    @property
    def sv_diagnostic_file(self) -> Path:
        return self.data_dir / "diagnostics" / "svs.tsv"

    def load(self, snv: bool = True, sv: bool = True) -> "CurationMaps":
        """Read the SNV and/or SV files; the others stay empty."""
        if snv:
            self.snv_validation = load_validation_map(self.snv_validation_file, None)
            self.snv_diagnostic = load_diagnostic_map(self.snv_diagnostic_file, None)
        if sv:
            self.sv_validation = load_validation_map(self.sv_validation_file, None)
            self.sv_diagnostic = load_diagnostic_map(self.sv_diagnostic_file, None)
        return self


def sv_variant_key(row: Dict[str, Any]) -> str:
    """The ``chr:start-end:dup|del`` key SV validations and diagnostics use.

    Built from the uncurated locus, as the SVs tab and the SV dialog build it,
    so a diagnostic saved here is read back by every table.
    """
    locus = row.get("_original_locus") or row.get("chr:start-end", "")
    return f"{locus}:{infer_sv_type(row)}"


def apply_curated_coordinates(
    row: Dict[str, Any],
    validation_map: Dict,
    sv_key: str,
    sample_id: str,
) -> None:
    """Apply curated coordinate bounds from validation data to a WCX row.

    If a 'present' validation with curated start/end exists, updates
    chr:start-end and Variant with curated values, sets IsCurated=True,
    and adds a tooltip showing original vs curated coordinates.
    """
    row["IsCurated"] = False
    row["_curated_tooltip"] = ""
    map_key = (sv_key, sample_id)
    if map_key not in validation_map:
        return
    validations = validation_map[map_key]
    present_with_curated = [
        v for v in validations if v[0] == "present" and v[3] != "1" and (v[4] or v[5])
    ]
    if not present_with_curated:
        return
    present_with_curated.sort(key=lambda v: v[6], reverse=True)
    most_recent = present_with_curated[0]
    curated_start = most_recent[4]
    curated_end = most_recent[5]

    original_locus = row.get("_original_locus", row.get("chr:start-end", ""))
    parts = original_locus.split(":")
    if len(parts) != 2:
        return
    chrom = parts[0]
    range_parts = parts[1].split("-")
    if len(range_parts) != 2:
        return
    orig_start = range_parts[0]
    orig_end = range_parts[1]
    new_start = curated_start if curated_start else orig_start
    new_end = curated_end if curated_end else orig_end
    row["chr:start-end"] = f"{chrom}:{new_start}-{new_end}"
    row["Variant"] = row["chr:start-end"]
    row["IsCurated"] = True
    row["_curated_tooltip"] = (
        f"Original: {original_locus}\nCurated: {row['chr:start-end']}"
    )
    # Recompute svlen from curated boundaries
    try:
        row["svlen"] = int(new_end) - int(new_start)
    except (ValueError, TypeError):
        pass


def refresh_curation(rows: List[Dict[str, Any]], maps: CurationMaps) -> None:
    """Re-read validations and diagnostics into rows already built."""
    for row in rows:
        s_id = row.get("sample", "")
        if row.get("_source_type") == "wisecondorx":
            orig = row.get("_original_locus", row.get("chr:start-end", ""))
            sv_vk = sv_variant_key(row)
            row["chr:start-end"] = orig
            row["Variant"] = orig
            add_validation_status_to_row(row, maps.sv_validation, sv_vk, s_id)
            add_diagnostic_status_to_row(row, maps.sv_diagnostic, sv_vk, s_id)
            apply_curated_coordinates(row, maps.sv_validation, sv_vk, s_id)
        else:
            v_key = row.get("Variant", "")
            add_validation_status_to_row(row, maps.snv_validation, v_key, s_id)
            add_diagnostic_status_to_row(row, maps.snv_diagnostic, v_key, s_id)


# ---------------------------------------------------------------------------
# Result rows
# ---------------------------------------------------------------------------


def _person_fields(
    row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, str]],
    sample_to_family: Dict[str, str],
) -> None:
    sample_id = row.get("sample", "")
    ped_info = pedigree_data.get(sample_id, {})
    row["FID"] = ped_info.get("FID", sample_to_family.get(sample_id, ""))
    row["Phenotype"] = ped_info.get("Phenotype", "")


def wombat_row_passes(row: Dict[str, Any], comps: List[Dict[str, Any]]) -> bool:
    """Check if a wombat row passes at least one component's filters.

    A row passes if ANY wombat component would keep it.
    """
    return any(
        consequences_selected(row, wc["impacts"]["value"])
        and not (wc["exclude_lcr"]["value"] and in_lcr(row))
        and not (wc["exclude_gnomad"]["value"] and gnomad_filtered(row))
        and not (wc["exclude_gnomad_wgs"]["value"] and in_gnomad_wgs(row))
        and status_selected(row.get("Validation"), wc["validations"]["value"])
        for wc in comps
    )


def wcx_row_passes(row: Dict[str, Any], comps: List[Dict[str, Any]]) -> bool:
    """Check if a WCX row passes at least one component's validation filter."""
    return any(
        status_selected(
            row.get("Validation"), wc.get("validations", {}).get("value", [])
        )
        for wc in comps
    )


def build_wombat_rows(
    frame: pl.DataFrame,
    comps: List[Dict[str, Any]],
    *,
    cohort_name: str,
    pedigree_data: Dict[str, Dict[str, str]],
    sample_to_family: Dict[str, str],
    maps: CurationMaps,
    gene_scorer: Any,
    unknown_consequences: set,
    unknown_clinvar: set,
) -> List[Dict[str, Any]]:
    """Wombat result rows with badges and curation, filtered by the components."""
    rows = frame.to_dicts()
    for row in rows:
        row["_source_type"] = "wombat"
        row["IsCurated"] = False
        row["_curated_tooltip"] = ""
        chrom = row.get("#CHROM", "")
        pos = row.get("POS", "")
        ref = row.get("REF", "")
        alt = row.get("ALT", "")
        sample_id = row.get("sample", "")
        variant_key = f"{chrom}:{pos}:{ref}:{alt}"
        row["Variant"] = variant_key
        row["_cohort_name"] = cohort_name
        _person_fields(row, pedigree_data, sample_to_family)

        add_vep_badges(row, unknown_consequences, unknown_clinvar)

        # Gene badges from VEP_SYMBOL (score-sorted, capped)
        row["GeneBadges"] = build_gene_badges(
            split_genes(row.get("VEP_SYMBOL", "")), gene_scorer
        )
        # Synthetic "gene" column for unified display
        row["gene"] = ",".join(badge_labels(row["GeneBadges"]))

        # Gene badges from VEP_Gene (ENSG IDs)
        gene_str = row.get("VEP_Gene", "")
        if gene_str:
            genes = [g.strip() for g in str(gene_str).split(",") if g.strip()]
            row["VEP_Gene_badges"] = build_gene_badges(
                genes, gene_scorer, cap=None, sort=False
            )
        else:
            row["VEP_Gene_badges"] = []

        add_validation_status_to_row(row, maps.snv_validation, variant_key, sample_id)
        add_diagnostic_status_to_row(row, maps.snv_diagnostic, variant_key, sample_id)

        add_score_badges(row)

    return [r for r in rows if wombat_row_passes(r, comps)]


def build_wcx_rows(
    frame: pl.DataFrame,
    comps: List[Dict[str, Any]],
    *,
    cohort_name: str,
    pedigree_data: Dict[str, Dict[str, str]],
    sample_to_family: Dict[str, str],
    maps: CurationMaps,
    gene_scorer: Any,
) -> List[Dict[str, Any]]:
    """WisecondorX result rows with badges, curation and curated bounds."""
    rows = frame.to_dicts()
    for row in rows:
        row["_source_type"] = "wisecondorx"
        sample_id = row.get("sample", "")
        row["_original_locus"] = row.get("chr:start-end", "")
        row["Variant"] = row.get("chr:start-end", "")
        row["_cohort_name"] = cohort_name
        # Compute svlen from start/end
        try:
            row["svlen"] = int(row.get("end", 0)) - int(row.get("start", 0))
        except (ValueError, TypeError):
            row["svlen"] = None
        _person_fields(row, pedigree_data, sample_to_family)

        # WisecondorX tooltip (ratio + zscore)
        tooltip_parts = []
        if row.get("ratio", "") not in (None, ""):
            tooltip_parts.append(f"ratio: {row.get('ratio')}")
        if row.get("zscore", "") not in (None, ""):
            tooltip_parts.append(f"zscore: {row.get('zscore')}")
        row["_wisecondorx_tooltip"] = ", ".join(tooltip_parts)

        # No VEP badges for SVS rows
        row["ConsequenceBadges"] = []
        row["ClinVarBadges"] = []

        # Gene badges from genic_symbol (score-sorted, capped)
        genic_str = row.get("genic_symbol", "")
        if genic_str and str(genic_str).strip():
            symbols = [s.strip() for s in str(genic_str).split(",") if s.strip()]
            row["GeneBadges"] = build_gene_badges(symbols, gene_scorer)
        else:
            row["GeneBadges"] = []
        # Synthetic "gene" column for unified display
        row["gene"] = ",".join(badge_labels(row["GeneBadges"]))

        row["VEP_Gene_badges"] = []

        # SVS validation: variant_key is chr:start-end:type
        sv_key = sv_variant_key(row)
        add_validation_status_to_row(row, maps.sv_validation, sv_key, sample_id)
        add_diagnostic_status_to_row(row, maps.sv_diagnostic, sv_key, sample_id)
        apply_curated_coordinates(row, maps.sv_validation, sv_key, sample_id)

    return [r for r in rows if wcx_row_passes(r, comps)]


# ---------------------------------------------------------------------------
# Columns, presets, sorting, export
# ---------------------------------------------------------------------------


def result_columns(
    wombat: Optional[pl.DataFrame], wcx: Optional[pl.DataFrame]
) -> List[str]:
    """Every column of the results, before grouping."""
    all_col_set: set[str] = set()
    if wombat is not None:
        all_col_set.update(wombat.columns)
    if wcx is not None:
        all_col_set.update(wcx.columns)
    # Remove internal/coordinate columns kept for filtering
    all_col_set -= {"_source_type", "chr", "start", "end"}
    all_col_set.add("gene")
    if wombat is not None:
        all_col_set.add("MOI")
    if wcx is not None:
        all_col_set.add("svlen")
    all_columns = list(all_col_set)
    if "Variant" not in all_columns:
        all_columns.insert(0, "Variant")
    for ensure_col in ["FID", "Phenotype", "Validation", "Diagnostic"]:
        if ensure_col not in all_columns:
            all_columns.append(ensure_col)
    return all_columns


def preset_named(presets: List[Dict[str, Any]], name: str) -> Dict[str, Any]:
    """The preset of that name, else the first one."""
    return next((p for p in presets if p["name"] == name), presets[0])


def initial_preset(
    components: List[Dict[str, Any]],
    source_files: List[Dict[str, Any]],
    presets: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Auto-select preset based on component mix.

    WCX only: SV View. Mixed: Mix View. Wombat only: by the first Wombat
    component's config name.
    """
    has_wombat = any(c["type"] == "wombat" for c in components)
    has_wcx = any(c["type"] == "wisecondorx" for c in components)
    if has_wcx and not has_wombat:
        return preset_named(presets, "SV View")
    if has_wcx and has_wombat:
        return preset_named(presets, "Mix View")
    first_wombat = next((c for c in components if c["type"] == "wombat"), None)
    wombat_config = ""
    if first_wombat:
        sf = find_source(source_files, first_wombat["source_key"])
        if sf:
            wombat_config = sf.get("wombat_config", "")
    return select_preset_for_config(wombat_config, presets)


def validation_file_columns(
    preset_name: str, all_columns: List[str], presets: List[Dict[str, Any]]
) -> List[str]:
    """Columns of a "Save as validation file" TSV.

    The SV and Mix views keep their own columns plus what the validation
    file page needs of an SV; every other view saves the Standard View.
    """
    sv_view = preset_name in ("SV View", "Mix View")
    target = preset_named(presets, preset_name if sv_view else "Standard View")
    save_cols = [
        col
        for col in target.get("columns", [])
        if col != "Validation" and col in all_columns and is_export_column(col)
    ]
    # Auto-include essential SV fields for validation file page
    if sv_view:
        for sv_col in ["ratio", "zscore", "type"]:
            if sv_col in all_columns and sv_col not in save_cols:
                save_cols.append(sv_col)
    return save_cols


def sort_rows(
    rows: List[Dict[str, Any]], sort_field: str, sort_type: str, desc: bool
) -> None:
    """Sort rows in place as the table would, so a restored sort starts right.

    ``DataTable`` does not reorder its initial rows to match
    ``initial_sorting``.
    """
    key: Callable[[Dict[str, Any]], Tuple[bool, Any]]
    if sort_type == "genomic":

        def key(r):
            return (r.get(sort_field) is None, genomic_sort_key(r.get(sort_field, "")))

    elif sort_type == "numerical":

        def key(r):
            v = r.get(sort_field)
            if v is None:
                return (True, 0.0)
            try:
                return (False, float(v))
            except (ValueError, TypeError):
                return (True, 0.0)

    else:

        def key(r):
            return (r.get(sort_field) is None, r.get(sort_field, ""))

    rows.sort(key=key, reverse=desc)
