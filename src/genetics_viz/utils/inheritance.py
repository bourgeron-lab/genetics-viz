"""Inheritance: what a trio's genotypes show, and whether it fits the gene.

Two layers:

- :func:`infer_inheritance` reads the transmission of one variant from the
  child's and parents' genotypes (de novo, paternal, ...). It pre-fills the
  validation dialog.
- :func:`assess_moi` compares that observation with the gene's known modes of
  inheritance from GeneTrek (AD, AR, XLD, XLR, ...) and gives a verdict for the
  MOI column: consistent, possible or inconsistent. :func:`annotate_moi` does
  it for a whole table, which is where recessive second hits are found.

The verdict is a reading aid, not a filter. Its limits are stated in each
tooltip: second hits are looked for only among the rows of the same table, and
GeneTrek's combined ``inheritance`` holds one mode per gene, which is why the
modes of every source are pooled (:func:`gene_modes`).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

from genetics_viz.utils.gene_scoring import split_genes
from genetics_viz.utils.genetrek import GenetrekSnapshot
from genetics_viz.utils.pedigree_labels import normalize_pheno, normalize_sex

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "inheritance.yaml"
with open(_CONFIG_PATH) as _f:
    INHERITANCE_CONFIG: Dict[str, Any] = yaml.safe_load(_f) or {}

_MODE_MAP: Dict[str, Optional[str]] = INHERITANCE_CONFIG.get("modes") or {}
MODE_LABELS: Dict[str, str] = INHERITANCE_CONFIG.get("mode_labels") or {}
STATUSES: Dict[str, Dict[str, Any]] = INHERITANCE_CONFIG.get("statuses") or {}
_PAR: List[Tuple[int, int]] = [
    (int(a), int(b)) for a, b in INHERITANCE_CONFIG.get("par_regions_grch38") or []
]
_warned: set = set()

CONSISTENT = "consistent"
POSSIBLE = "possible"
INCONSISTENT = "inconsistent"


# ---------------------------------------------------------------------------
# Genotypes and transmission
# ---------------------------------------------------------------------------

# Genotype values understood by infer_inheritance (unphased form).
_GT_CARRIES = {"0/1", "1/0", "1/1"}
_GT_ABSENT = {"0/0"}
_GT_MISSING = {"./.", ".", ""}


def normalize_gt(gt: Optional[str]) -> str:
    """One genotype string: phasing dropped, ``""`` for a disagreeing list.

    Aggregated rows join distinct values with ``,`` (``"0/1,1/1"`` when two
    annotation rows of one variant disagree); such a value means nothing.
    """
    if gt is None:
        return ""
    text = str(gt).strip().replace("|", "/")
    if "," in text:
        parts = {p.strip() for p in text.split(",") if p.strip()}
        if len(parts) != 1:
            return ""
        text = parts.pop()
    return text


def infer_inheritance(sample_gt: str, father_gt: str, mother_gt: str) -> str:
    """Infer inheritance mode from VCF-style genotype strings.

    Args:
        sample_gt: Child/sample genotype (e.g. "0/1", "1/1", "0/0")
        father_gt: Father's genotype
        mother_gt: Mother's genotype

    Returns:
        One of the inheritance dropdown values:
        "unknown", "de novo", "paternal", "maternal",
        "not paternal", "not maternal", "either", "homozygous"
    """
    s = normalize_gt(sample_gt)
    f = normalize_gt(father_gt)
    m = normalize_gt(mother_gt)

    f_carries = f in _GT_CARRIES
    f_absent = f in _GT_ABSENT
    f_missing = f in _GT_MISSING

    m_carries = m in _GT_CARRIES
    m_absent = m in _GT_ABSENT
    m_missing = m in _GT_MISSING

    # --- Hom-alt sample (1/1) ---
    if s == "1/1":
        # Both parents carry → homozygous (autosomal recessive)
        if f_carries and m_carries:
            return "homozygous"
        return "unknown"

    # --- Het sample (0/1 or 1/0) ---
    if s in ("0/1", "1/0"):
        # Both parents known
        if f_carries and m_carries:
            return "either"
        if f_carries and m_absent:
            return "paternal"
        if f_absent and m_carries:
            return "maternal"
        if f_absent and m_absent:
            return "de novo"

        # One parent known, other missing
        if f_carries and m_missing:
            return "paternal"
        if m_carries and f_missing:
            return "maternal"
        if f_absent and m_missing:
            return "not paternal"
        if m_absent and f_missing:
            return "not maternal"

        return "unknown"

    # --- All other cases (0/0, ./., missing, etc.) ---
    return "unknown"


def _alt_count(gt: str) -> Optional[int]:
    """Alt alleles in a normalised genotype; ``None`` when not called."""
    alleles = gt.split("/") if gt else []
    if not alleles or any(a in (".", "") for a in alleles):
        return None
    return sum(a != "0" for a in alleles)


def _ploidy(gt: str) -> int:
    return len(gt.split("/")) if gt else 0


def in_par(chrom: str, pos: Optional[int]) -> bool:
    return pos is not None and any(a <= pos <= b for a, b in _PAR)


def _chrom(chrom: Optional[str]) -> str:
    c = str(chrom or "").strip()
    return c[3:] if c.lower().startswith("chr") else c


@dataclass(frozen=True)
class Observation:
    """What one row shows for one sample."""

    #: "het", "hom", "hemi" or "unknown".
    zygosity: str
    #: "de novo", "paternal", "maternal", "both", "either" or "unknown".
    origin: str
    #: Pedigree phenotype ("1" unaffected, "2" affected, "-9") of the parent
    #: the variant came from, when exactly one is known.
    parent_pheno: str = "-9"
    on_x: bool = False
    #: Normalised pedigree sex: "1" male, "2" female, "-9" unknown.
    sex: str = "-9"

    def describe(self) -> str:
        zyg = {
            "het": "heterozygous",
            "hom": "homozygous",
            "hemi": "hemizygous",
        }.get(self.zygosity, "genotype unclear")
        origin = {
            "both": "from both parents",
            "either": "both parents carry it",
            "paternal": "paternal",
            "maternal": "maternal",
            "de novo": "de novo",
        }.get(self.origin, "origin unknown")
        return f"{zyg}, {origin}"


def observe(
    *,
    chrom: Optional[str],
    pos: Any,
    sample_gt: Optional[str],
    father_gt: Optional[str] = None,
    mother_gt: Optional[str] = None,
    sex: Optional[str] = None,
    father_pheno: Optional[str] = None,
    mother_pheno: Optional[str] = None,
    de_novo: bool = False,
) -> Observation:
    """Read zygosity and parental origin from one row's genotypes.

    ``de_novo`` marks a row from a de novo caller, which carries no genotypes:
    heterozygous (or hemizygous) by construction.
    """
    c = _chrom(chrom)
    try:
        position = int(pos) if pos is not None else None
    except (TypeError, ValueError):
        position = None
    sex_code = normalize_sex(str(sex) if sex is not None else None)
    on_x = c.upper() == "X"
    haploid_x = on_x and sex_code == "1" and not in_par(c, position)

    if de_novo:
        return Observation(
            "hemi" if haploid_x else "het", "de novo", on_x=on_x, sex=sex_code
        )

    s, f, m = normalize_gt(sample_gt), normalize_gt(father_gt), normalize_gt(mother_gt)
    n_alt = _alt_count(s)
    if not n_alt:
        return Observation("unknown", "unknown", on_x=on_x, sex=sex_code)
    if haploid_x or (_ploidy(s) == 1 and on_x):
        zygosity = "hemi"
    elif n_alt >= 2:
        zygosity = "hom"
    else:
        zygosity = "het"

    f_alt, m_alt = _alt_count(f), _alt_count(m)
    f_carries, m_carries = bool(f_alt), bool(m_alt)
    f_absent, m_absent = f_alt == 0, m_alt == 0

    if zygosity == "hemi":
        # A male's X comes from his mother.
        origin = "maternal" if m_carries else "de novo" if m_absent else "unknown"
    elif zygosity == "hom":
        origin = "both" if f_carries and m_carries else "unknown"
    elif f_carries and m_carries:
        origin = "either"
    elif f_carries and not m_carries:
        origin = "paternal"
    elif m_carries and not f_carries:
        origin = "maternal"
    elif f_absent and m_absent:
        origin = "de novo"
    else:
        origin = "unknown"

    parent_pheno = "-9"
    if origin == "paternal":
        parent_pheno = normalize_pheno(father_pheno)
    elif origin == "maternal":
        parent_pheno = normalize_pheno(mother_pheno)
    return Observation(zygosity, origin, parent_pheno, on_x, sex_code)


# ---------------------------------------------------------------------------
# Gene modes
# ---------------------------------------------------------------------------


def mode_code(raw: str) -> Optional[str]:
    """GeneTrek's raw mode string -> a code, via ``inheritance.yaml``."""
    if raw in _MODE_MAP:
        return _MODE_MAP[raw]
    if raw not in _warned:
        _warned.add(raw)
        logger.warning(
            "Unknown GeneTrek inheritance value %r, add it to inheritance.yaml", raw
        )
    return None


def gene_modes(
    snapshot: GenetrekSnapshot, symbol: str
) -> Tuple[List[str], List[Tuple[str, str]]]:
    """The pooled mode codes of a gene, and each source's raw value."""
    info = snapshot.info.get(symbol)
    sources = list(info.inheritance) if info else []
    codes: List[str] = []
    for _, raw in sources:
        code = mode_code(raw)
        if code and code not in codes:
            codes.append(code)
    return codes, sources


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SecondHit:
    """Another heterozygous variant of the same gene in the same sample."""

    variant: str
    origin: str


@dataclass(frozen=True)
class Verdict:
    status: str
    mode: str
    reason: str

    @property
    def rank(self) -> int:
        return int(STATUSES[self.status]["rank"])


def _assess_mode(
    mode: str, obs: Observation, second_hits: Sequence[SecondHit]
) -> Optional[Verdict]:
    zyg, origin = obs.zygosity, obs.origin
    if mode == "MT" or zyg == "unknown":
        return None

    if mode in ("AD", "XLD"):
        if origin in ("paternal", "maternal") and obs.parent_pheno == "1":
            return Verdict(
                POSSIBLE,
                mode,
                f"inherited from an unaffected {_parent(origin)}: reduced penetrance?",
            )
        if origin == "de novo":
            return Verdict(CONSISTENT, mode, "de novo, as expected for a dominant gene")
        return Verdict(CONSISTENT, mode, "one copy is enough for a dominant gene")

    if mode == "XLR" or (mode == "XL" and obs.on_x):
        if zyg == "hemi":
            return Verdict(CONSISTENT, mode, "hemizygous male")
        if zyg == "hom":
            return Verdict(CONSISTENT, mode, "homozygous")
        if obs.sex == "1":
            return Verdict(CONSISTENT, mode, "male, het call on chrX")
        recessive = _recessive_het(mode, second_hits, origin)
        if recessive is not None and recessive.status == CONSISTENT:
            return recessive
        who = "female" if obs.sex == "2" else "sex unknown"
        return Verdict(
            POSSIBLE, mode, f"heterozygous, {who}: carrier, or skewed X-inactivation?"
        )

    if mode == "AR":
        if zyg in ("hom", "hemi"):
            if origin == "both":
                return Verdict(CONSISTENT, mode, "homozygous, both parents carriers")
            return Verdict(
                POSSIBLE,
                mode,
                "homozygous, parents not both seen carrying it: deletion, UPD or missing parent?",
            )
        return _recessive_het(mode, second_hits, origin)

    return None


def _parent(origin: str) -> str:
    return "father" if origin == "paternal" else "mother"


def _recessive_het(mode: str, second_hits: Sequence[SecondHit], origin: str) -> Verdict:
    if not second_hits:
        return Verdict(
            INCONSISTENT, mode, "a single heterozygous variant in this table"
        )
    parental = {"paternal", "maternal"}
    others = {h.origin for h in second_hits}
    if origin in parental and (others & parental) - {origin}:
        return Verdict(
            CONSISTENT, mode, "compound heterozygous, in trans (one from each parent)"
        )
    if origin in parental and others == {origin}:
        return Verdict(
            INCONSISTENT, mode, f"two heterozygous variants, both {origin}: in cis"
        )
    return Verdict(POSSIBLE, mode, "two or more heterozygous variants, phase unknown")


def assess_moi(
    modes: Iterable[str], obs: Observation, second_hits: Sequence[SecondHit] = ()
) -> Optional[Verdict]:
    """The best verdict over the gene's modes, or ``None`` with nothing to say."""
    verdicts = [v for m in modes if (v := _assess_mode(m, obs, second_hits))]
    return max(verdicts, key=lambda v: v.rank) if verdicts else None


# ---------------------------------------------------------------------------
# A table at a time
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RowContext:
    """How to read a table's rows."""

    gene_column: str = "VEP_SYMBOL"
    sample_column: str = "sample"
    chrom_column: str = "#CHROM"
    pos_column: str = "POS"
    #: Rows from a de novo caller: no genotype columns.
    de_novo: bool = False


def family_pedigree(
    samples: Iterable[Any],
) -> Callable[[str], Dict[str, Optional[str]]]:
    """A pedigree lookup over ``models.Sample`` objects, for :func:`annotate_moi`."""
    by_id = {s.sample_id: s for s in samples}

    def lookup(sample_id: str) -> Dict[str, Optional[str]]:
        s = by_id.get(sample_id)
        if s is None:
            return {}
        father = by_id.get(s.father_id) if s.father_id else None
        mother = by_id.get(s.mother_id) if s.mother_id else None
        return {
            "sex": s.sex,
            "father_pheno": father.phenotype if father else None,
            "mother_pheno": mother.phenotype if mother else None,
        }

    return lookup


def _variant_label(row: Dict[str, Any], ctx: RowContext) -> str:
    ref, alt = row.get("REF"), row.get("ALT")
    label = f"{row.get(ctx.chrom_column, '')}:{row.get(ctx.pos_column, '')}"
    return f"{label} {ref}>{alt}" if ref and alt else label


def annotate_moi(
    rows: List[Dict[str, Any]],
    snapshot: GenetrekSnapshot,
    *,
    pedigree: Callable[[str], Dict[str, Optional[str]]],
    ctx: RowContext = RowContext(),
) -> None:
    """Add ``MoiBadges`` and ``_moi_rank`` to each row, in place.

    ``pedigree(sample)`` returns the sample's ``sex`` and its parents'
    phenotypes as ``father_pheno`` / ``mother_pheno`` (values may be ``None``).
    """
    observed: List[Tuple[Observation, List[str], str]] = []
    hets: Dict[Tuple[str, str], List[Tuple[int, SecondHit]]] = {}
    for i, row in enumerate(rows):
        sample = str(row.get(ctx.sample_column) or "")
        ped = pedigree(sample) if sample else {}
        obs = observe(
            chrom=row.get(ctx.chrom_column),
            pos=row.get(ctx.pos_column),
            sample_gt=row.get("sample_gt"),
            father_gt=row.get("father_gt"),
            mother_gt=row.get("mother_gt"),
            sex=row.get("sex") or ped.get("sex"),
            father_pheno=ped.get("father_pheno"),
            mother_pheno=ped.get("mother_pheno"),
            de_novo=ctx.de_novo,
        )
        genes: List[str] = []
        for ident in split_genes(row.get(ctx.gene_column)):
            symbol = snapshot.resolve(ident)
            if symbol and symbol not in genes:
                genes.append(symbol)
        observed.append((obs, genes, sample))
        if obs.zygosity == "het":
            hit = SecondHit(_variant_label(row, ctx), obs.origin)
            for symbol in genes:
                hets.setdefault((sample, symbol), []).append((i, hit))

    for i, (row, (obs, genes, sample)) in enumerate(zip(rows, observed)):
        badges: List[Dict[str, Any]] = []
        best = 0
        for symbol in genes:
            modes, sources = gene_modes(snapshot, symbol)
            if not modes:
                continue
            others = [h for j, h in hets.get((sample, symbol), []) if j != i]
            verdict = assess_moi(modes, obs, others)
            if verdict is None:
                continue
            best = max(best, verdict.rank)
            badges.append(
                _badge(symbol, len(genes) > 1, verdict, obs, modes, sources, others)
            )
        row["MoiBadges"] = badges
        row["_moi_rank"] = best
        row["MOI"] = ", ".join(b["label"] for b in badges)


def _badge(
    symbol: str,
    prefix: bool,
    verdict: Verdict,
    obs: Observation,
    modes: List[str],
    sources: List[Tuple[str, str]],
    others: List[SecondHit],
) -> Dict[str, Any]:
    style = STATUSES[verdict.status]
    label = f"{verdict.mode} {style['label']}"
    if prefix:
        label = f"{symbol} {label}"
    lines = [
        f"{symbol}: {verdict.status} with {MODE_LABELS.get(verdict.mode, verdict.mode)}",
        f"Observed: {obs.describe()}",
        f"Reading: {verdict.reason}",
        "Modes in GeneTrek: " + ", ".join(modes),
    ]
    lines += [f"  {name}: {raw}" for name, raw in sources]
    if others:
        listed = "; ".join(f"{h.variant} ({h.origin})" for h in others[:5])
        more = f" and {len(others) - 5} more" if len(others) > 5 else ""
        lines.append(f"Other het variants of {symbol} in this table: {listed}{more}")
    if verdict.mode == "AR" or verdict.mode.startswith("XL"):
        lines.append("Second hits are looked for only in the rows of this table.")
    return {
        "label": label,
        "color": style["color"],
        "textColor": style["text_color"],
        "tooltip": "\n".join(lines),
    }
