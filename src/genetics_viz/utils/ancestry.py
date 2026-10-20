"""Ancestry / polygenic-score file discovery, reference panel loading and colours.

The ``ancestry-pgs`` pipeline writes one ``ancestry/`` subdirectory per family,
holding files named ``<family_id>.apgs_b<bundle>_<filters>.<kind>.tsv`` - for
example ``88171200005.apgs_b1.0.0_dp10gq20.pcs.tsv``. The principal components
in those files are projections onto a *versioned* reference panel, so the panel
that a family must be plotted against is the one named by its own filename tag,
not a globally configured version. Only the root of the references tree is
configured (``ancestry_reference_dir``); the ``v<bundle>/`` subdirectory below
it is derived per family.
"""

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Collection, Dict, List, Optional, Tuple

import polars as pl
import yaml

from ..config_model import get_ancestry_reference_dir
from .export import PERSON_COLUMNS, people_frame
from .sharding import get_family_path
from .tsv import read_tsv_or_none

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "ancestry.yaml"

#: Sentinels used for missing values across the ancestry TSVs. The PGS z-score
#: file writes the literal string ``nan``, which polars would otherwise type as
#: a string and make the whole column unsortable.
NULL_VALUES = [".", "", "nan", "NaN", "NA"]

#: Identifier columns of the ancestry-pgs tables, always read as strings. The
#: cohort tables append ``family_id``, which is numeric for most families of
#: some cohorts and not for others.
ID_COLUMNS = ("IID", "FID", "family_id", "father", "mother")

#: Principal components kept from the reference panel. It carries PC1-PC20, but
#: only the first four are plotted and the panel has ~6.7k rows.
REFERENCE_PCS = ["PC1", "PC2", "PC3", "PC4"]

_FALLBACK_COLOR = "#94a3b8"


def _load_config() -> Dict[str, Any]:
    """Load the ancestry display configuration from YAML."""
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


_CONFIG: Dict[str, Any] = _load_config()

REGION_COLORS: Dict[str, str] = _CONFIG.get("region_colors", {})


def get_region_color(region: Optional[str]) -> str:
    """Return the display colour for a reference panel region."""
    if not region:
        return REGION_COLORS.get("unknown", _FALLBACK_COLOR)
    for key, color in REGION_COLORS.items():
        if key.lower() == region.lower():
            return color
    return REGION_COLORS.get("unknown", _FALLBACK_COLOR)


def get_pgs_zscore_thresholds() -> List[Dict[str, Any]]:
    """Return the ``color_scale`` threshold list for PGS z-score columns."""
    return list(_CONFIG.get("pgs_zscore_thresholds", []))


def reload_ancestry_config() -> None:
    """Reload the ancestry display configuration from YAML."""
    global _CONFIG, REGION_COLORS
    _CONFIG = _load_config()
    REGION_COLORS = _CONFIG.get("region_colors", {})


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

#: ``<entity_id>.apgs_b<bundle>_<filters>.<kind>.tsv``. The bundle version
#: contains dots, so it needs a lazy multi-segment capture rather than the
#: single dot-free segment used by the wombat/SV patterns.
_FILENAME_TMPL = (
    r"{eid}\.apgs_b(?P<bundle>.+?)_(?P<filters>[^.]+)\.(?P<kind>{kinds})\.tsv$"
)

#: Result kinds of the ancestry step, in its ``ancestry/`` directory. ``Q`` is
#: the admixture proportions (``IID``, ``snp_coverage``, one column per
#: reference component); only the exports read it.
_KINDS = ("pcs", "ancestry", "pgs_zscore", "Q")

#: Result kinds of the ROH step, in its own ``roh/`` directory. Same naming:
#: the step labels its outputs with the ancestry panel they were called on.
#: ``roh`` and ``roh_genes`` cannot be confused, the pattern is anchored on
#: ``.tsv``.
ROH_KINDS = ("froh", "roh", "roh_genes", "kinship", "roh_windows")


def _filename_pattern(prefix: str, kinds: Tuple[str, ...]) -> "re.Pattern[str]":
    """Compile the result-file pattern for one entity and a set of kinds."""
    return re.compile(
        _FILENAME_TMPL.format(
            eid=re.escape(prefix), kinds="|".join(re.escape(k) for k in kinds)
        )
    )


#: A ``(bundle_version, filter_tag)`` pair identifying one set of result files.
Variant = Tuple[str, str]


@dataclass
class AncestryFiles:
    """The ancestry result files of one family, for a single bundle version."""

    directory: Path
    bundle_version: str
    filter_tag: str
    paths: Dict[str, Path] = field(default_factory=dict)
    #: Every ``(bundle_version, filter_tag)`` pair present in the directory,
    #: best first, so a caller can offer the choice. Always includes the
    #: selected one.
    all_variants: List[Variant] = field(default_factory=list)

    @property
    def variant(self) -> Variant:
        """The selected ``(bundle_version, filter_tag)`` pair."""
        return (self.bundle_version, self.filter_tag)

    @property
    def other_variants(self) -> List[Variant]:
        """The variants present but not selected."""
        return [v for v in self.all_variants if v != self.variant]

    @property
    def pcs_path(self) -> Optional[Path]:
        return self.paths.get("pcs")

    @property
    def ancestry_path(self) -> Optional[Path]:
        return self.paths.get("ancestry")

    @property
    def pgs_zscore_path(self) -> Optional[Path]:
        return self.paths.get("pgs_zscore")

    def path(self, kind: str) -> Optional[Path]:
        """The file of ``kind`` for the selected variant, if present."""
        return self.paths.get(kind)

    @property
    def file_tag(self) -> str:
        """The filename tag, e.g. ``apgs_b1.0.0_dp10gq20``."""
        return f"apgs_b{self.bundle_version}_{self.filter_tag}"

    @property
    def label(self) -> str:
        """Human-readable provenance, e.g. ``bundle 1.0.0 / dp10gq20``."""
        return f"bundle {self.bundle_version} / {self.filter_tag}"


_VERSION_PART_RE = re.compile(r"^(\d*)(.*)$")


def _version_key(version: str) -> Tuple[Tuple[int, int, str], ...]:
    """Sort key for a dotted version string.

    Each part sorts on its leading number, then a release before anything with
    a suffix, so ``1.0.0`` ranks above ``1.0.0rc1`` - joining all the digits of
    a part would read ``0rc1`` as 1 and put the candidate above the release.
    Unparseable parts sort as 0.
    """
    parts = []
    for part in version.split("."):
        number, suffix = _VERSION_PART_RE.match(part).groups()
        parts.append((int(number) if number else 0, 0 if suffix else 1, suffix))
    return tuple(parts)


_TAG_RE = re.compile(r"^dp(?P<dp>\d+)gq(?P<gq>\d+)$")


def _tag_strictness(tag: str) -> Tuple[int, int, int, str]:
    """Sort key ranking a filter tag by how strictly it filtered genotypes.

    Tags look like ``dp10gq20`` - minimum depth 10, minimum GQ 20. Sorted
    descending, the strictest variant comes first, so ``dp10gq20`` is preferred
    over ``dp6gq15``. A plain string sort would pick ``dp6gq15`` instead,
    because "6" > "1" at the third character.

    Unrecognised tags sort below every recognised one and fall back to
    comparing the tag itself, keeping the order stable.
    """
    match = _TAG_RE.match(tag)
    if match is None:
        return (0, 0, 0, tag)
    return (1, int(match.group("dp")), int(match.group("gq")), tag)


def _variant_key(
    variant: Variant,
) -> Tuple[Tuple[Tuple[int, int, str], ...], Tuple[int, int, int, str]]:
    """Sort key for a variant: newest bundle first, then strictest filtering."""
    bundle, tag = variant
    return (_version_key(bundle), _tag_strictness(tag))


def get_ancestry_dir(data_dir: Path, family_id: str) -> Path:
    """Return the ``ancestry/`` directory of a family (may not exist)."""
    return get_family_path(data_dir, family_id) / "ancestry"


def _scan_variants(
    directory: Path, prefix: str, kinds: Tuple[str, ...] = _KINDS
) -> Dict[Variant, Dict[str, Path]]:
    """Map every ``(bundle, tag)`` variant in ``directory`` to its files."""
    pattern = _filename_pattern(prefix, kinds)
    variants: Dict[Variant, Dict[str, Path]] = {}
    for tsv_file in directory.glob("*.tsv"):
        match = pattern.match(tsv_file.name)
        if not match:
            continue
        key = (match.group("bundle"), match.group("filters"))
        variants.setdefault(key, {})[match.group("kind")] = tsv_file
    return variants


def find_ancestry_files_in(
    directory: Path,
    prefix: str,
    variant: Optional[Variant] = None,
    kinds: Tuple[str, ...] = _KINDS,
) -> Optional[AncestryFiles]:
    """Locate ancestry result files named ``<prefix>.apgs_b<bundle>_<tag>.*``.

    Used for both layouts: a family's ``ancestry/`` directory, where the prefix
    is the family ID, and a cohort's, where it is the cohort name. ``kinds``
    selects which step's files are wanted - :data:`ROH_KINDS` for the ``roh/``
    directories, which share the naming.

    Without an explicit ``variant``, the newest bundle and then the strictest
    filter tag wins - see :func:`_tag_strictness`. Every variant found is
    reported in :attr:`AncestryFiles.all_variants` so a caller can offer the
    choice. Returns ``None`` when the directory is missing, holds no recognised
    file, or does not hold the requested variant.
    """
    if not directory.is_dir():
        return None

    variants = _scan_variants(directory, prefix, kinds)
    if not variants:
        return None

    ordered = sorted(variants, key=_variant_key, reverse=True)
    if variant is None:
        selected = ordered[0]
    elif variant in variants:
        selected = variant
    else:
        return None

    bundle, tag = selected
    return AncestryFiles(
        directory=directory,
        bundle_version=bundle,
        filter_tag=tag,
        paths=variants[selected],
        all_variants=ordered,
    )


def find_ancestry_files(
    data_dir: Path, family_id: str, variant: Optional[Variant] = None
) -> Optional[AncestryFiles]:
    """Locate a family's ancestry result files."""
    return find_ancestry_files_in(
        get_ancestry_dir(data_dir, family_id), family_id, variant
    )


def get_cohort_ancestry_dir(cohort: Any) -> Path:
    """Return the ``ancestry/`` directory of a cohort (may not exist).

    Takes a :class:`~genetics_viz.models.Cohort`; typed loosely to keep this
    module free of a models import.
    """
    return cohort.path / "ancestry"


def find_cohort_ancestry_files(
    cohort: Any, variant: Optional[Variant] = None
) -> Optional[AncestryFiles]:
    """Locate a cohort's ancestry result files.

    Cohort-level files carry an extra ``family_id`` column that the per-family
    ones do not.
    """
    return find_ancestry_files_in(get_cohort_ancestry_dir(cohort), cohort.name, variant)


def probe_kind_in(directory: Path, prefix: str, kind: str) -> bool:
    """Return True when ``directory`` holds a ``<prefix>`` result file of ``kind``."""
    if not directory.is_dir():
        return False
    pattern = _filename_pattern(prefix, (kind,))
    return any(pattern.match(f.name) for f in directory.glob("*.tsv"))


def _probe(data_dir: Path, family_id: str, kind: str) -> bool:
    """Return True when a file of ``kind`` exists for this family."""
    return probe_kind_in(get_ancestry_dir(data_dir, family_id), family_id, kind)


def probe_ancestry_data(data_dir: Path, family_id: str) -> bool:
    """Check whether principal components exist for this family."""
    return _probe(data_dir, family_id, "pcs")


def probe_pgs_data(data_dir: Path, family_id: str) -> bool:
    """Check whether PGS z-scores exist for this family."""
    return _probe(data_dir, family_id, "pgs_zscore")


# ---------------------------------------------------------------------------
# Reference panel
# ---------------------------------------------------------------------------

#: ``(reference_root, bundle_version)`` -> panel frame, or None when absent.
#: The panel is ~6.7k rows and identical for every family, so it is read once
#: per bundle rather than on every page render.
_reference_cache: Dict[Tuple[str, str], Optional[pl.DataFrame]] = {}


def get_reference_pcs_path(bundle_version: str) -> Optional[Path]:
    """Return the reference PCs file for a bundle version, or ``None``.

    ``None`` means the ``ancestry_reference_dir`` config key is unset; the
    returned path is not guaranteed to exist.
    """
    root = get_ancestry_reference_dir()
    if root is None:
        return None
    return root / f"v{bundle_version}" / "bundle" / "labels" / "reference_pcs.tsv.gz"


def load_reference_pcs(bundle_version: str) -> Optional[pl.DataFrame]:
    """Load the reference panel PCs for a bundle version.

    Returns a frame of ``IID, PC1-PC4, population, region``, or ``None`` when
    the reference directory is not configured, the file is missing, or it
    cannot be parsed. A missing panel is a degraded-but-usable state - the
    caller still plots the family's own samples - so it is not an error.
    """
    path = get_reference_pcs_path(bundle_version)
    if path is None:
        return None

    cache_key = (str(path.parent), bundle_version)
    if cache_key in _reference_cache:
        return _reference_cache[cache_key]

    frame: Optional[pl.DataFrame] = None
    if not path.is_file():
        logger.warning(
            "Reference panel not found for bundle %s: %s", bundle_version, path
        )
    else:
        try:
            df = read_tsv_or_none(
                path, null_values=NULL_VALUES, string_columns=ID_COLUMNS
            )
            if df is None or len(df) == 0:
                logger.warning("Reference panel is empty: %s", path)
            else:
                wanted = ["IID", *REFERENCE_PCS, "population", "region"]
                missing = [c for c in wanted if c not in df.columns]
                if missing:
                    logger.warning(
                        "Reference panel %s is missing columns: %s", path, missing
                    )
                else:
                    frame = df.select(wanted)
        except Exception:
            logger.warning("Failed to read reference panel %s", path, exc_info=True)

    _reference_cache[cache_key] = frame
    return frame


def clear_reference_cache() -> None:
    """Drop the cached reference panels (used when config or data changes)."""
    _reference_cache.clear()


# ---------------------------------------------------------------------------
# Exports, shared by the family and cohort tabs
# ---------------------------------------------------------------------------

#: Columns of the pipeline tables that are keys, not data.
_KEY_COLUMNS = ("IID", "FID", "family_id")


def _read_export(path: Optional[Path]) -> Optional[pl.DataFrame]:
    if path is None:
        return None
    df = read_tsv_or_none(
        path,
        infer_schema_length=10000,
        null_values=NULL_VALUES,
        string_columns=ID_COLUMNS,
    )
    if df is None or "IID" not in df.columns:
        return None
    return df


def _attach(
    base: pl.DataFrame,
    frame: Optional[pl.DataFrame],
    how: str = "left",
    rename: Optional[Dict[str, str]] = None,
) -> Tuple[pl.DataFrame, List[str]]:
    """Join ``frame``'s data columns onto ``base``; return them with the result.

    On ``(family_id, IID)`` when the frame carries ``family_id`` (the cohort
    tables: a sample may sit in two families), on ``IID`` otherwise.
    """
    if frame is None:
        return (base.clear() if how == "inner" else base), []
    data = [c for c in frame.columns if c not in _KEY_COLUMNS and c not in base.columns]
    if "family_id" in frame.columns:
        right = frame.select(["family_id", "IID", *data]).rename({"family_id": "FID"})
        on = ["FID", "IID"]
    else:
        right = frame.select(["IID", *data])
        on = ["IID"]
    right = right.unique(subset=on, keep="first", maintain_order=True)
    if rename:
        right = right.rename({k: v for k, v in rename.items() if k in data})
        data = [rename.get(c, c) for c in data]
    return base.join(right, on=on, how=how, maintain_order="left"), data


def ancestry_export(files: AncestryFiles, people: List[Dict[str, Any]]) -> pl.DataFrame:
    """The people with principal components: labels, PCs, admixture proportions.

    Blocking I/O. The admixture components are prefixed ``admix_``; a family
    without a ``Q.tsv`` gets none.
    """
    base = people_frame(people)
    out, pcs = _attach(base, _read_export(files.pcs_path), how="inner")
    out, labels = _attach(out, _read_export(files.ancestry_path))
    q = _read_export(files.path("Q"))
    components = (
        []
        if q is None
        else [c for c in q.columns if c not in _KEY_COLUMNS and c != "snp_coverage"]
    )
    out, admixture = _attach(out, q, rename={c: f"admix_{c}" for c in components})
    return out.select([*PERSON_COLUMNS, *labels, *pcs, *admixture])


def pgs_export(
    files: AncestryFiles,
    people: List[Dict[str, Any]],
    regions: Optional[Collection[str]] = None,
) -> pl.DataFrame:
    """The people with z-scores: region, then every trait.

    ``regions`` keeps only those predicted regions (``"—"`` for none), as the
    cohort tab's region selector does.
    """
    base = people_frame(people)
    out, traits = _attach(base, _read_export(files.pgs_zscore_path), how="inner")
    labels = _read_export(files.ancestry_path)
    if labels is not None and "region" in labels.columns:
        keys = [c for c in ("IID", "family_id") if c in labels.columns]
        out, _ = _attach(out, labels.select([*keys, "region"]))
    else:
        out = out.with_columns(pl.lit(None, dtype=pl.String).alias("region"))
    if regions is not None:
        out = out.filter(pl.col("region").fill_null("—").is_in(list(regions)))
    return out.select([*PERSON_COLUMNS, "region", *traits])
