"""Normalised coverage and the CNV control panel, as written by ghfc-ngs.

The ``cnv_coverage`` step writes, per sample,
``samples/.../<ID>/coverage/<ID>.cov_b<bin>q<mapq>.{bw,json}``: every bin of a
fixed grid holds ``2 x depth / autosomal median``, so a normal autosome reads
about 2. The ``cnv_panel`` step of a dedicated controls cohort writes a
versioned, self-contained directory ``cohorts/<controls>/cnv_panel/v<version>/``
with the per-bin median and unscaled MAD of the controls for ``autosomes``,
``XX`` (chrX) and ``XY`` (chrX, chrY), and ``panel.json``. A cohort fitted
against a panel also gets
``cohorts/<NAME>/coverage/<NAME>.<tag>.panel_v<version>.fit.tsv``.

The division happens here, at view time: corrected copy number = expected
ploidy x sample / panel median. The per-sample files never depend on a panel,
so a new panel version needs no rewrite of any sample.

Every missing piece is a degraded state, never an error. Without a panel, or
with one of another tag, the sample's own normalised value is drawn (about the
copy number on autosomes), unmasked, with the reason said.
"""

from __future__ import annotations

import json
import logging
import math
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import yaml

from genetics_viz.utils.pedigree_labels import normalize_sex
from genetics_viz.utils.sharding import get_sample_path
from genetics_viz.utils.tsv import read_tsv_or_none

logger = logging.getLogger(__name__)


_CONFIG_PATH = Path(__file__).parent.parent / "config" / "cnv.yaml"


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH) as f:
        return yaml.safe_load(f) or {}


#: ``config/cnv.yaml``, shared with ``sv_view`` and ``sv_reads``.
CNV_CONFIG: Dict[str, Any] = _load_config()
CN_CONFIG: Dict[str, Any] = CNV_CONFIG.get("copy_number", {})

COVERAGE_SUBDIR = "coverage"
TAG_RE = re.compile(r"^cov_b(\d+)q(\d+)$")
_SAMPLE_FILE_RE = re.compile(r"^(?P<sample>.+)\.(?P<tag>cov_b\d+q\d+)\.bw$")

AUTOSOMES = frozenset(f"chr{i}" for i in range(1, 23))
#: Karyotypes the panel has a group for. Anything else the pipeline spells
#: out (X0, XXY, ambiguous) gets the raw value on the sex chromosomes.
SEX_GROUPS = ("XX", "XY")


# ---------------------------------------------------------------------------
# Sample files
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SampleCoverage:
    sample_id: str
    tag: str
    bigwig: Path
    json_path: Path


def find_sample_coverage(data_dir: Path, sample_id: str) -> List[SampleCoverage]:
    """Return the sample's coverage bigWigs, one per tag, sorted by tag."""
    cov_dir = get_sample_path(data_dir, sample_id) / COVERAGE_SUBDIR
    if not cov_dir.is_dir():
        return []
    found = []
    for path in cov_dir.iterdir():
        m = _SAMPLE_FILE_RE.match(path.name)
        if m and m.group("sample") == sample_id:
            found.append(
                SampleCoverage(
                    sample_id=sample_id,
                    tag=m.group("tag"),
                    bigwig=path,
                    json_path=path.with_suffix(".json"),
                )
            )
    return sorted(found, key=lambda c: c.tag)


def pick_coverage(
    coverages: Sequence[SampleCoverage], tag: Optional[str]
) -> Optional[SampleCoverage]:
    """The coverage of the panel's tag when there is one, else the first."""
    if not coverages:
        return None
    if tag:
        for cov in coverages:
            if cov.tag == tag:
                return cov
    return coverages[0]


def autosomal_median_depth(coverage: SampleCoverage) -> Optional[float]:
    """The sample's ``autosomal_median_depth`` from its coverage JSON."""
    try:
        value = json.loads(coverage.json_path.read_text()).get("autosomal_median_depth")
        depth = float(value)
    except (OSError, ValueError, TypeError, AttributeError):
        return None
    return depth if depth > 0 else None


# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Panel:
    directory: Path
    version: str
    tag: str
    bin: int
    max_fit_mad: Optional[float]
    source_cohort: str
    n_controls: Dict[str, int]
    #: ``chrom -> [(start0, end)]``, half-open 0-based, from 1-based inclusive.
    par: Dict[str, List[Tuple[int, int]]]
    expected_ploidy: Dict[str, Any]

    def bigwig(self, group: str, stat: str) -> Path:
        return self.directory / f"panel.{self.tag}.{group}.{stat}.bw"

    @property
    def label(self) -> str:
        return (
            f"{self.source_cohort} v{self.version}"
            if self.source_cohort
            else f"v{self.version}"
        )


#: ``(path, mtime, size)`` of ``panel.json`` -> parsed panel (``None`` = invalid).
_panel_cache: Dict[Tuple[str, float, int], Optional[Panel]] = {}


def panel_dir_from_params(raw: Any, data_dirs: Sequence[Path]) -> Optional[Path]:
    """Map the params file's ``cnv_panel_dir`` onto a local directory.

    The value is a compute-cluster path (``/pasteur/helix/.../cohorts/
    GHFC-NT/cnv_panel/v1``) that does not resolve here. Its tail from the last
    ``cohorts/`` component is resolved under each data directory in turn, the
    cohort's own first. The cluster path itself is never opened.
    """
    if not raw or not isinstance(raw, str):
        return None
    parts = Path(raw.strip()).parts
    try:
        idx = len(parts) - 1 - parts[::-1].index("cohorts")
    except ValueError:
        return None
    tail = Path(*parts[idx:])
    for data_dir in data_dirs:
        candidate = Path(data_dir) / tail
        if (candidate / "panel.json").is_file():
            return candidate
    return None


def _parse_par(raw: Any) -> Dict[str, List[Tuple[int, int]]]:
    par: Dict[str, List[Tuple[int, int]]] = {}
    if not isinstance(raw, dict):
        return par
    for chrom, regions in raw.items():
        if not isinstance(regions, dict):
            continue  # "build", "coordinates"
        spans = []
        for bounds in regions.values():
            if isinstance(bounds, (list, tuple)) and len(bounds) == 2:
                spans.append((int(bounds[0]) - 1, int(bounds[1])))
        par[str(chrom)] = sorted(spans)
    return par


def load_panel(panel_dir: Path) -> Optional[Panel]:
    """Read a panel version directory, cached on ``panel.json``'s stat."""
    path = Path(panel_dir) / "panel.json"
    try:
        stat = path.stat()
    except OSError:
        return None
    key = (str(path), stat.st_mtime, stat.st_size)
    if key in _panel_cache:
        return _panel_cache[key]
    for old in [k for k in _panel_cache if k[0] == key[0]]:
        del _panel_cache[old]

    panel: Optional[Panel] = None
    try:
        data = json.loads(path.read_text())
        tag = str(data["tag"])
        if not TAG_RE.match(tag):
            raise ValueError(f"unexpected tag {tag!r}")
        max_fit = data.get("max_fit_mad")
        panel = Panel(
            directory=Path(panel_dir),
            version=str(data["panel_version"]),
            tag=tag,
            bin=int(data["bin"]),
            max_fit_mad=float(max_fit) if max_fit is not None else None,
            source_cohort=str(data.get("source_cohort") or ""),
            n_controls={
                str(k): int(v) for k, v in (data.get("n_controls") or {}).items()
            },
            par=_parse_par(data.get("par")),
            expected_ploidy=dict(data["expected_ploidy"]),
        )
        if not panel.bigwig("autosomes", "median").is_file():
            raise ValueError("no autosomal median bigWig")
    except (OSError, ValueError, KeyError, TypeError) as e:
        logger.warning("Ignoring CNV panel %s: %s", panel_dir, e)
        panel = None
    _panel_cache[key] = panel
    return panel


def _own_panel_dir(cohort_dir: Path, params: Dict[str, Any]) -> Optional[Path]:
    """The panel a controls cohort built itself: its ``cnv_panel_version``, else
    the newest complete version under ``cnv_panel/``."""
    root = Path(cohort_dir) / "cnv_panel"
    version = params.get("cnv_panel_version")
    if version:
        candidate = root / f"v{version}"
        if (candidate / "panel.json").is_file():
            return candidate
    if not root.is_dir():
        return None
    versions = [d for d in root.iterdir() if (d / "panel.json").is_file()]
    return max(versions, key=lambda d: d.stat().st_mtime) if versions else None


def resolve_cohort_panel(
    params: Dict[str, Any],
    data_dirs: Sequence[Path],
    cohort_dir: Optional[Path] = None,
) -> Optional[Panel]:
    """The panel a cohort is divided by.

    The one its params' ``cnv_panel_dir`` names. Failing that, a controls
    cohort uses the panel it built itself, whose params name no
    ``cnv_panel_dir``.
    """
    panel_dir = panel_dir_from_params(params.get("cnv_panel_dir"), data_dirs)
    if panel_dir is None and cohort_dir is not None:
        panel_dir = _own_panel_dir(cohort_dir, params)
    return load_panel(panel_dir) if panel_dir else None


def clear_cnv_cache() -> None:
    """Drop the cached panels and fit tables."""
    _panel_cache.clear()
    _fit_cache.clear()


# ---------------------------------------------------------------------------
# Cohort fit
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FitRow:
    family_id: str
    sample: str
    pedigree_sex: str
    coverage_sex: str
    fit_mad: Optional[float]
    chrX_ratio: Optional[float]
    chrY_ratio: Optional[float]


_fit_cache: Dict[Tuple[str, float, int], Dict[Tuple[str, str], FitRow]] = {}


def fit_file(cohort_dir: Path, cohort_name: str, panel: Panel) -> Path:
    return (
        Path(cohort_dir)
        / COVERAGE_SUBDIR
        / f"{cohort_name}.{panel.tag}.panel_v{panel.version}.fit.tsv"
    )


def _float_or_none(value: Any) -> Optional[float]:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def read_cohort_fit(
    cohort_dir: Path, cohort_name: str, panel: Optional[Panel]
) -> Dict[Tuple[str, str], FitRow]:
    """The cohort's fit against ``panel``, keyed by ``(family_id, sample)``.

    Only the table of this panel's tag and version is read: a fit against
    another panel says nothing about this one.
    """
    if panel is None:
        return {}
    path = fit_file(cohort_dir, cohort_name, panel)
    try:
        stat = path.stat()
    except OSError:
        return {}
    key = (str(path), stat.st_mtime, stat.st_size)
    cached = _fit_cache.get(key)
    if cached is not None:
        return cached
    for old in [k for k in _fit_cache if k[0] == key[0]]:
        del _fit_cache[old]

    rows: Dict[Tuple[str, str], FitRow] = {}
    try:
        df = read_tsv_or_none(path, string_columns=("family_id", "sample"))
    except Exception as e:  # a malformed table is a degraded state
        logger.warning("Could not read CNV fit table %s: %s", path, e)
        df = None
    if df is not None:
        for r in df.iter_rows(named=True):
            fid, sample = str(r.get("family_id") or ""), str(r.get("sample") or "")
            rows[(fid, sample)] = FitRow(
                family_id=fid,
                sample=sample,
                pedigree_sex=str(r.get("pedigree_sex") or ""),
                coverage_sex=str(r.get("coverage_sex") or ""),
                fit_mad=_float_or_none(r.get("fit_mad")),
                chrX_ratio=_float_or_none(r.get("chrX_ratio")),
                chrY_ratio=_float_or_none(r.get("chrY_ratio")),
            )
    _fit_cache[key] = rows
    return rows


def karyotype(fit: Optional[FitRow], pedigree_sex: Optional[str]) -> Optional[str]:
    """``"XX"``, ``"XY"`` or ``None`` when no panel group fits the sample.

    The coverage sex from the fit table wins: it is what the pipeline measured.
    One it spells otherwise (``X0``, ``XXY``, ``ambiguous``) has no group. With
    no fit row, the pedigree sex stands in.
    """
    if fit is not None and fit.coverage_sex:
        return fit.coverage_sex if fit.coverage_sex in SEX_GROUPS else None
    code = normalize_sex(pedigree_sex)
    return {"1": "XY", "2": "XX"}.get(code)


# ---------------------------------------------------------------------------
# Copy number
# ---------------------------------------------------------------------------


@dataclass
class CNTrack:
    sample_id: str
    #: Display-bin centres (bp) and values; ``None`` where there is no data.
    positions: List[int] = field(default_factory=list)
    values: List[Optional[float]] = field(default_factory=list)
    masked: List[bool] = field(default_factory=list)
    #: ``"corrected"`` (divided by the panel) or ``"raw"`` (own normalisation).
    mode: str = "raw"
    reason: Optional[str] = None
    karyotype: Optional[str] = None
    fit_mad: Optional[float] = None
    bin_bp: int = 0


def _panel_group(chrom: str, kary: Optional[str]) -> Optional[str]:
    if chrom in AUTOSOMES:
        return "autosomes"
    if chrom == "chrX" and kary in SEX_GROUPS:
        return kary
    if chrom == "chrY" and kary == "XY":
        return "XY"
    return None


def expected_ploidy(
    panel: Panel, group: str, chrom: str, centres: np.ndarray
) -> np.ndarray:
    """Expected copy number per bin, PAR-aware, from ``panel.json``."""
    table = panel.expected_ploidy
    if group == "autosomes":
        return np.full(centres.shape, float(table.get("autosomes", 2)))
    entry = (table.get(group) or {}).get(chrom)
    if isinstance(entry, (int, float)):
        return np.full(centres.shape, float(entry))
    if not isinstance(entry, dict):
        return np.full(centres.shape, np.nan)
    out = np.full(centres.shape, float(entry.get("non_par", np.nan)))
    if "par" in entry:
        for s, e in panel.par.get(chrom, []):
            out[(centres >= s) & (centres < e)] = float(entry["par"])
    else:
        # chrY of XY has no PAR ploidy: its PARs are N-masked and read on chrX.
        for s, e in panel.par.get(chrom, []):
            out[(centres >= s) & (centres < e)] = np.nan
    return out


def _resolve_chrom(chroms: Dict[str, int], chrom: str) -> Optional[str]:
    if chrom in chroms:
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    return alt if alt in chroms else None


#: Per-call memo of bigWig reads, so the panel is read once for a whole trio.
ReadMemo = Dict[Tuple[str, str, int, int, int, bool], np.ndarray]


def _read(
    path: Path,
    chrom: str,
    start: int,
    end: int,
    bins: int,
    exact: bool,
    memo: Optional[ReadMemo] = None,
) -> np.ndarray:
    key = (str(path), chrom, start, end, bins, exact)
    if memo is not None and key in memo:
        return memo[key]
    out = _read_bigwig(path, chrom, start, end, bins, exact)
    if memo is not None:
        memo[key] = out
    return out


def _read_bigwig(
    path: Path, chrom: str, start: int, end: int, bins: int, exact: bool
) -> np.ndarray:
    import pybigtools

    with pybigtools.open(str(path)) as bw:
        if chrom not in bw.chroms():
            return np.full(bins, np.nan)
        return np.asarray(
            bw.values(
                chrom, start, end, bins=bins, summary="mean", exact=exact, fillna=None
            ),
            dtype=float,
        )


def _aggregate(
    values: np.ndarray, masked: np.ndarray, n_points: int
) -> Tuple[np.ndarray, np.ndarray, int]:
    """Median-aggregate consecutive bins to at most ``n_points``.

    A display bin is masked when fewer than half of its bins are usable; its
    value is then the median of all of them, for a grey point in the right
    place rather than a gap.
    """
    n = len(values)
    k = max(1, math.ceil(n / max(1, n_points)))
    if k == 1:
        return values, masked, 1
    m = math.ceil(n / k)
    pad = m * k - n
    v = np.concatenate([values, np.full(pad, np.nan)]).reshape(m, k)
    mk = np.concatenate([masked, np.ones(pad, dtype=bool)]).reshape(m, k)
    usable = np.where(mk, np.nan, v)
    n_usable = np.sum(~np.isnan(usable), axis=1)
    with warnings.catch_warnings():
        # nanmedian of an all-NaN slice warns; an empty bin is ordinary here.
        warnings.simplefilter("ignore", RuntimeWarning)
        good = np.nanmedian(usable, axis=1)
        anyv = np.nanmedian(v, axis=1)
    out_masked = n_usable < (k / 2)
    out = np.where(out_masked, anyv, good)
    return out, out_masked, k


def copy_number(
    coverage: SampleCoverage,
    panel: Optional[Panel],
    kary: Optional[str],
    chrom: str,
    start: int,
    end: int,
    *,
    fit_mad: Optional[float] = None,
    cfg: Optional[Dict[str, Any]] = None,
    memo: Optional[ReadMemo] = None,
) -> CNTrack:
    """Copy number of one sample over ``chrom:start-end``, ready to draw.

    Blocking (bigWig reads over the network mount): call it off-thread. Pass
    the same ``memo`` for every sample of one view, so the panel is read once.
    """
    import pybigtools

    cfg = cfg or CN_CONFIG
    track = CNTrack(sample_id=coverage.sample_id, karyotype=kary, fit_mad=fit_mad)
    m = TAG_RE.match(coverage.tag)
    bin_bp = int(m.group(1)) if m else 250

    with pybigtools.open(str(coverage.bigwig)) as bw:
        chroms = dict(bw.chroms())
    bw_chrom = _resolve_chrom(chroms, chrom)
    if bw_chrom is None:
        track.reason = f"{chrom} is not in the coverage file"
        return track
    chrom_len = chroms[bw_chrom]
    # Whole bins of the pipeline's grid only, so the reads line up with it;
    # the truncated last bin of a contig is left out.
    s0 = max(0, (int(start) // bin_bp) * bin_bp)
    e0 = min((chrom_len // bin_bp) * bin_bp, -(-int(end) // bin_bp) * bin_bp)
    if e0 <= s0:
        return track
    n_bins = (e0 - s0) // bin_bp
    max_exact = int(cfg.get("max_exact_bins", 200_000))
    n_points = int(cfg.get("display_points", 1500))
    exact = n_bins <= max_exact
    nb = n_bins if exact else n_points
    width = (e0 - s0) / nb
    centres = s0 + (np.arange(nb) + 0.5) * width

    sample = _read(coverage.bigwig, bw_chrom, s0, e0, nb, exact)  # own file, no memo

    group = _panel_group(bw_chrom, kary)
    reason = None
    if panel is None:
        reason = "no control panel: own normalisation, not divided"
    elif panel.tag != coverage.tag:
        reason = f"panel is {panel.tag}, sample is {coverage.tag}: not divided"
    elif group is None:
        reason = f"no panel group for karyotype {kary or 'unknown'} on {bw_chrom}: not divided"

    if reason is not None:
        cn = sample
        masked = np.isnan(sample)
        track.mode = "raw"
    else:
        assert panel is not None and group is not None
        median = _read(panel.bigwig(group, "median"), bw_chrom, s0, e0, nb, exact, memo)
        mad_path = panel.bigwig(group, "mad")
        mad = (
            _read(mad_path, bw_chrom, s0, e0, nb, exact, memo)
            if mad_path.is_file()
            else np.zeros(nb)
        )
        expected = expected_ploidy(panel, group, bw_chrom, centres)
        with np.errstate(all="ignore"):
            cn = expected * sample / median
            rel_mad = mad / median
        min_frac = float(cfg.get("min_median_fraction", 0.5))
        max_rel = float(cfg.get("max_relative_mad", 0.35))
        masked = (
            np.isnan(cn)
            | np.isnan(expected)
            | ~(median > 0)
            | (median < min_frac * expected)
            | (rel_mad > max_rel)
        )
        # Where the panel says nothing, show the sample's own value, greyed.
        cn = np.where(np.isfinite(cn), cn, sample)
        track.mode = "corrected"

    values, out_masked, k = _aggregate(cn, masked, n_points)
    if k > 1:
        centres = s0 + (np.arange(len(values)) + 0.5) * width * k
    y_max = float(cfg.get("y_max", 5))
    track.positions = [int(c) for c in centres[: len(values)]]
    track.values = [
        None if not np.isfinite(v) else round(float(min(v, y_max)), 3) for v in values
    ]
    track.masked = [bool(x) for x in out_masked]
    track.reason = reason
    track.bin_bp = int(round(width * k))
    return track


# ---------------------------------------------------------------------------
# One view
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CNRequest:
    sample_id: str
    family_id: str
    pedigree_sex: Optional[str]


def cn_window(
    start: int, end: int, cfg: Optional[Dict[str, Any]] = None
) -> Tuple[int, int]:
    """The SV span times the dialog's expansion, at least ``min_window_bp``."""
    cfg = cfg or CN_CONFIG
    expansion = float(CNV_CONFIG.get("read_view", {}).get("expansion", 1.6))
    span = max(1, int(end) - int(start))
    width = max(int(span * expansion), int(cfg.get("min_window_bp", 20000)))
    centre = (int(start) + int(end)) // 2
    lo = max(0, centre - width // 2)
    return lo, lo + width


def load_cn_view(
    data_dir: Path,
    data_dirs: Sequence[Path],
    cohort_dir: Optional[Path],
    cohort_name: str,
    params: Dict[str, Any],
    requests: Sequence[CNRequest],
    chrom: str,
    start: int,
    end: int,
) -> Tuple[Optional[Panel], List[Optional[CNTrack]]]:
    """Copy-number tracks of several samples over one window.

    Blocking: call it off-thread. A sample with no coverage file gets
    ``None``. The panel is that of the cohort's params, read once.
    """
    panel = resolve_cohort_panel(params, data_dirs, cohort_dir)
    fits = read_cohort_fit(cohort_dir, cohort_name, panel) if cohort_dir else {}
    by_sample = {k[1]: v for k, v in fits.items()}
    memo: ReadMemo = {}
    tracks: List[Optional[CNTrack]] = []
    for req in requests:
        cov = pick_coverage(
            find_sample_coverage(data_dir, req.sample_id), panel.tag if panel else None
        )
        if cov is None:
            tracks.append(None)
            continue
        fit = fits.get((req.family_id, req.sample_id)) or by_sample.get(req.sample_id)
        tracks.append(
            copy_number(
                cov,
                panel,
                karyotype(fit, req.pedigree_sex),
                chrom,
                start,
                end,
                fit_mad=fit.fit_mad if fit else None,
                memo=memo,
            )
        )
    return panel, tracks
