"""The ghfc-ngs ``denovo`` step: discovery, reading, evidence on Wombat rows, filters.

The step annotates every candidate of the Wombat de novo table (its
``#wombat_config``, ``de_novo_mutations`` by default) with region masks, a
pileup of every member's CRAM, parent of origin (Unfazed) and flags, and counts
each child's *strict* candidates against the callable genome. It writes:

- family: ``families/…/<FID>/denovo/<FID>.<label>.{dnm,qc,inputs}.tsv``;
- cohort: ``cohorts/<NAME>/denovo/<NAME>.<label>.{dnm,qc,inputs}.tsv``, the
  family tables concatenated, plus ``recurrence_families`` and the
  ``RECURRENT`` flag (which can turn a strict candidate flagged), and the QC
  verdict per child (``qc``: ok, high, low, no_callable) with its range and the
  trio's sampleqc flags.

The pipeline decides every flag and status; the app words them from
``config/denovo.yaml`` and filters on them. The family table joins the Wombat
rows 1:1 on (chrom, pos, ref, alt, child): the step reads that table and writes
one row per candidate. Where the cohort table is as new as the family's, its
rows are used, so the family page shows ``RECURRENT`` too.

The QC table's ``n_strict_snv``/``n_strict_indel``/``ti_fraction`` are the
family pass's, written before ``RECURRENT``; the cohort tab recounts them from
the candidates so every figure follows one strict set.
"""

from __future__ import annotations

import csv
import logging
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import polars as pl
import yaml

from .sharding import get_family_path
from .tsv import find_labels, read_tsv_or_none
from .validation_status import IN_PHASE_MNV, consensus

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "denovo.yaml"

DENOVO_SUBDIR = "denovo"
#: The Wombat config the step reads unless its inputs say otherwise.
DEFAULT_WOMBAT_CONFIG = "de_novo_mutations"
#: The app's own flag: a Wombat row the step has no evidence for.
NOT_ASSESSED = "NOT_ASSESSED"
AUTOSOMES = frozenset(f"chr{i}" for i in range(1, 23))

_GREY = "#9e9e9e"

_STRING_COLUMNS = (
    "family_id", "child", "father", "mother", "child_sex", "chrom", "ref", "alt",
    "var_type", "inheritance_model", "genes", "consequences", "wombat_child_gt",
    "sibs_called", "sibs_with_alt", "origin", "origin_status", "flags", "strict",
    "callable_status", "sampleqc_flags", "qc_reference", "qc",
)  # fmt: skip


def _load_config() -> Dict[str, Any]:
    try:
        with open(_CONFIG_PATH) as fh:
            return yaml.safe_load(fh) or {}
    except Exception:
        logger.warning("Cannot read %s", _CONFIG_PATH, exc_info=True)
        return {}


_CONFIG = _load_config()
STRICT_IGNORED: Set[str] = set(_CONFIG.get("strict_ignored") or [])
CATEGORIES: Dict[str, Dict[str, str]] = _CONFIG.get("categories") or {}
FLAGS: Dict[str, Dict[str, str]] = _CONFIG.get("flags") or {}
#: Flags in display order (the config's).
FLAG_ORDER: List[str] = list(FLAGS)
REFERENCES: List[Dict[str, Any]] = _CONFIG.get("references") or []
EXPECTED: Dict[str, Any] = _CONFIG.get("expected") or {}

#: Inheritances that say a validated call is not de novo.
INHERITED = frozenset({"paternal", "maternal", "either", "homozygous"})


# ---------------------------------------------------------------------------
# Wording
# ---------------------------------------------------------------------------


def flag_info(flag: str) -> Dict[str, str]:
    """Label, colour, tooltip and category of a flag; an unknown one raw, grey."""
    entry = FLAGS.get(flag)
    if entry is None:
        return {"label": flag, "color": _GREY, "tooltip": flag, "category": ""}
    category = entry.get("category", "")
    color = entry.get("color") or CATEGORIES.get(category, {}).get("color", _GREY)
    return {
        "label": entry.get("label", flag),
        "color": color,
        "tooltip": entry.get("tooltip", flag),
        "category": category,
    }


def flag_rank(flag: str) -> int:
    return FLAG_ORDER.index(flag) if flag in FLAG_ORDER else len(FLAG_ORDER)


def split_flags(value: Any) -> List[str]:
    """``"A,B"`` -> ``["A", "B"]``; ``.``, empty and null give no flag."""
    if value is None:
        return []
    return [f for f in str(value).split(",") if f and f != "."]


def flag_badges(flags: Iterable[str]) -> List[Dict[str, str]]:
    """``badge_list`` entries, in display order, each with its explanation."""
    badges = []
    for flag in sorted(set(flags), key=flag_rank):
        info = flag_info(flag)
        badges.append(
            {
                "label": info["label"],
                "color": info["color"],
                "tooltip": f"{flag}\n{info['tooltip']}",
            }
        )
    return badges


def _entry(section: str, key: Optional[str]) -> Dict[str, str]:
    return (_CONFIG.get(section) or {}).get(key or "", {}) or {}


def strict_badges(strict: Optional[str]) -> List[Dict[str, str]]:
    entry = _entry("strict", strict)
    if not entry:
        return []
    return [
        {
            "label": entry["label"],
            "color": entry["color"],
            "tooltip": entry.get("tooltip", ""),
        }
    ]


def origin_badges(
    origin: Optional[str], status: Optional[str], evidence: Any
) -> List[Dict[str, str]]:
    """The parent of origin, or nothing; the tooltip says how it was reached."""
    status_text = (_CONFIG.get("origin_status") or {}).get(status or "", status or "")
    entry = _entry("origin", origin)
    if not entry:
        return []
    tip = status_text
    if evidence not in (None, ""):
        tip += f" (Unfazed evidence: {evidence})"
    return [{"label": entry["label"], "color": entry["color"], "tooltip": tip}]


def origin_status_text(status: Optional[str]) -> str:
    return (_CONFIG.get("origin_status") or {}).get(status or "", status or "")


def qc_info(value: Optional[str]) -> Dict[str, str]:
    entry = _entry("qc", value)
    if not entry:
        return {"label": value or "–", "color": _GREY, "tooltip": value or ""}
    return {
        "label": entry["label"],
        "color": entry["color"],
        "tooltip": entry.get("tooltip", ""),
    }


def qc_reference_text(value: Optional[str]) -> str:
    return (_CONFIG.get("qc_reference") or {}).get(value or "", value or "")


# ---------------------------------------------------------------------------
# Discovery and reading
# ---------------------------------------------------------------------------


def get_family_denovo_dir(data_dir: Path, family_id: str) -> Path:
    return get_family_path(data_dir, family_id) / DENOVO_SUBDIR


def get_cohort_denovo_dir(cohort: Any) -> Path:
    return cohort.path / DENOVO_SUBDIR


def denovo_path(directory: Path, prefix: str, label: str, kind: str) -> Path:
    """``<directory>/<prefix>.<label>.<kind>.tsv``."""
    return directory / f"{prefix}.{label}.{kind}.tsv"


def find_denovo_labels(directory: Path, prefix: str) -> List[str]:
    """Labels with a candidates table in ``directory``, the newest written first."""
    return find_labels(directory, prefix, "dnm")


def probe_family_denovo(data_dir: Path, family_id: str) -> bool:
    return bool(
        find_denovo_labels(get_family_denovo_dir(data_dir, family_id), family_id)
    )


def probe_cohort_denovo(cohort: Any) -> bool:
    return bool(find_denovo_labels(get_cohort_denovo_dir(cohort), cohort.name))


def read_denovo_table(path: Optional[Path]) -> Optional[pl.DataFrame]:
    """One dnm or qc table: ``.`` as null, identifiers as strings.

    ``None`` when the path is ``None``, missing, or empty.
    """
    if path is None or not path.is_file():
        return None
    return read_tsv_or_none(
        path,
        null_values=[".", ""],
        infer_schema_length=10000,
        string_columns=_STRING_COLUMNS,
        quote_char=None,
    )


def read_inputs_meta(path: Path) -> Dict[str, str]:
    """The ``#key<TAB>value`` lines of an inputs snapshot (wombat_config, masks…)."""
    meta: Dict[str, str] = {}
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                parts = line[1:].rstrip("\n").split("\t", 1)
                if len(parts) == 2:
                    meta[parts[0]] = parts[1]
    except OSError:
        pass
    return meta


StatKey = Tuple[str, float, int]
_cache: Dict[StatKey, Optional[pl.DataFrame]] = {}


def _stat(path: Path) -> Optional[StatKey]:
    try:
        st = path.stat()
    except OSError:
        return None
    return (str(path), st.st_mtime, st.st_size)


def _read_cached(path: Path) -> Optional[pl.DataFrame]:
    key = _stat(path)
    if key is None:
        return None
    if key not in _cache:
        for old in [k for k in _cache if k[0] == key[0]]:
            del _cache[old]
        try:
            _cache[key] = read_denovo_table(path)
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            _cache[key] = None
    return _cache[key]


def clear_denovo_cache() -> None:
    """Drop every cached table (tests)."""
    _cache.clear()


@dataclass
class DenovoData:
    """The denovo tables of a family or a cohort, for one label."""

    label: str
    labels: List[str]
    dnm: pl.DataFrame
    qc: Optional[pl.DataFrame]
    meta: Dict[str, str] = field(default_factory=dict)
    #: ``cohort`` when the rows are the cohort table's (``RECURRENT`` known).
    level: str = "family"

    @property
    def wombat_config(self) -> str:
        return self.meta.get("wombat_config") or DEFAULT_WOMBAT_CONFIG

    @property
    def has_recurrence(self) -> bool:
        return "recurrence_families" in self.dnm.columns

    def rows(self) -> List[Dict[str, Any]]:
        return self.dnm.to_dicts()

    def qc_rows(self) -> List[Dict[str, Any]]:
        return self.qc.to_dicts() if self.qc is not None else []

    def qc_of(self, child: str) -> Optional[Dict[str, Any]]:
        return next((q for q in self.qc_rows() if q.get("child") == child), None)


def _empty_dnm() -> pl.DataFrame:
    return pl.DataFrame({c: [] for c in ("family_id", "child", "chrom", "pos")})


def load_cohort_denovo(
    cohort: Any, label: Optional[str] = None
) -> Optional[DenovoData]:
    """The cohort tables for ``label`` (default: the newest). Blocking, cached."""
    directory = get_cohort_denovo_dir(cohort)
    labels = find_denovo_labels(directory, cohort.name)
    if not labels:
        return None
    chosen = label if label in labels else labels[0]
    dnm = _read_cached(denovo_path(directory, cohort.name, chosen, "dnm"))
    return DenovoData(
        label=chosen,
        labels=labels,
        dnm=dnm if dnm is not None else _empty_dnm(),
        qc=_read_cached(denovo_path(directory, cohort.name, chosen, "qc")),
        meta=read_inputs_meta(denovo_path(directory, cohort.name, chosen, "inputs")),
        level="cohort",
    )


def load_family_denovo(
    data_dir: Path, family_id: str, cohort: Any = None, label: Optional[str] = None
) -> Optional[DenovoData]:
    """A family's denovo tables, from the cohort's when they are as new.

    The cohort table carries ``RECURRENT`` (and the strict status that
    follows) and the per-child QC verdict. It is used when it has this label
    and this family's rows and was written after the family table; otherwise
    the family's own tables. ``None`` when the family has none.
    """
    directory = get_family_denovo_dir(data_dir, family_id)
    labels = find_denovo_labels(directory, family_id)
    if not labels:
        return None
    chosen = label if label in labels else labels[0]
    fam_dnm_path = denovo_path(directory, family_id, chosen, "dnm")
    meta = read_inputs_meta(denovo_path(directory, family_id, chosen, "inputs"))
    if cohort is not None:
        cdir = get_cohort_denovo_dir(cohort)
        c_path = denovo_path(cdir, cohort.name, chosen, "dnm")
        c_stat, f_stat = _stat(c_path), _stat(fam_dnm_path)
        if c_stat is not None and (f_stat is None or c_stat[1] >= f_stat[1]):
            c_dnm = _read_cached(c_path)
            if c_dnm is not None and "family_id" in c_dnm.columns:
                mine = c_dnm.filter(pl.col("family_id") == family_id)
                if len(mine):
                    c_qc = _read_cached(denovo_path(cdir, cohort.name, chosen, "qc"))
                    if c_qc is not None and "family_id" in c_qc.columns:
                        c_qc = c_qc.filter(pl.col("family_id") == family_id)
                    return DenovoData(chosen, labels, mine, c_qc, meta, "cohort")
    dnm = _read_cached(fam_dnm_path)
    return DenovoData(
        label=chosen,
        labels=labels,
        dnm=dnm if dnm is not None else _empty_dnm(),
        qc=_read_cached(denovo_path(directory, family_id, chosen, "qc")),
        meta=meta,
        level="family",
    )


def denovo_source(
    wombat_files: List[Dict[str, Any]], data: Optional[DenovoData]
) -> Optional[Dict[str, Any]]:
    """The Wombat file the De novo sub-tab shows: the step's, else a de novo one.

    The step names its config in its inputs; without its tables, the
    ``de_novo_mutations`` config, else the first config matching ``de_?novo``.
    """
    from .wombat import DENOVO_CONFIG

    by_name = {wf["wombat_config"]: wf for wf in wombat_files}
    if data is not None and data.wombat_config in by_name:
        return by_name[data.wombat_config]
    if DEFAULT_WOMBAT_CONFIG in by_name:
        return by_name[DEFAULT_WOMBAT_CONFIG]
    return next(
        (wf for wf in wombat_files if DENOVO_CONFIG.search(wf["wombat_config"])), None
    )


def load_child_qc(
    data_dir: Path, family_id: str, cohort: Any = None
) -> List[Dict[str, Any]]:
    """The family's per-child QC rows, light: the QC tables only (a few KB).

    The verdict (``qc``) is in the cohort table only; the family table is the
    fallback, without it. For the family header, read at page build.
    """
    directory = get_family_denovo_dir(data_dir, family_id)
    labels = find_denovo_labels(directory, family_id)
    if not labels:
        return []
    if cohort is not None:
        c_qc = _read_cached(
            denovo_path(get_cohort_denovo_dir(cohort), cohort.name, labels[0], "qc")
        )
        if c_qc is not None and "family_id" in c_qc.columns:
            mine = c_qc.filter(pl.col("family_id") == family_id)
            if len(mine):
                return mine.to_dicts()
    qc = _read_cached(denovo_path(directory, family_id, labels[0], "qc"))
    return qc.to_dicts() if qc is not None else []


def qc_alerts(qc_rows: Iterable[Dict[str, Any]]) -> List[Tuple[str, str]]:
    """``(child, verdict)`` of the children whose de novo QC is not ``ok``."""
    return [
        (str(q.get("child")), str(q.get("qc")))
        for q in qc_rows
        if q.get("qc") and q.get("qc") != "ok"
    ]


# ---------------------------------------------------------------------------
# Evidence on the Wombat rows
# ---------------------------------------------------------------------------


def candidate_key(row: Dict[str, Any]) -> Tuple[str, str]:
    """``(chrom:pos:ref:alt, child)`` of a dnm row, the Wombat rows' key."""
    return (
        f"{row.get('chrom')}:{row.get('pos')}:{row.get('ref')}:{row.get('alt')}",
        str(row.get("child") or ""),
    )


#: The evidence columns the de novo sub-tab adds, in table order.
EVIDENCE_COLUMNS = [
    "dnm_strict", "dnm_flags", "dnm_child_reads", "dnm_child_vaf", "dnm_father_alt",
    "dnm_mother_alt", "dnm_sibs", "dnm_origin", "dnm_recurrence", "dnm_cluster_n",
    "dnm_nearby_calls", "dnm_binom_p", "dnm_strand_p",
]  # fmt: skip

#: badge_list columns: column -> (badges field, sort field)
EVIDENCE_BADGES = {
    "dnm_strict": ("DnmStrictBadges", None),
    "dnm_flags": ("DnmFlagBadges", "_dnm_n_flags"),
    "dnm_origin": ("DnmOriginBadges", None),
}


def _int(value: Any) -> Optional[int]:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _float(value: Any) -> Optional[float]:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if x == x else None


def set_evidence(row: Dict[str, Any], ev: Optional[Dict[str, Any]]) -> None:
    """The ``dnm_*`` fields of one Wombat row from its dnm row (or none)."""
    if ev is None:
        flags = [NOT_ASSESSED]
        row.update(
            dnm_strict="",
            dnm_flags=NOT_ASSESSED,
            DnmFlagBadges=flag_badges(flags),
            DnmStrictBadges=[],
            DnmOriginBadges=[],
            _dnm_n_flags=1,
            _dnm=None,
        )
        for col in EVIDENCE_COLUMNS:
            row.setdefault(col, None)
        return
    flags = split_flags(ev.get("flags"))
    alt, ref = _int(ev.get("child_alt")), _int(ev.get("child_ref"))
    sibs = ev.get("sibs_with_alt") or ev.get("sibs_called") or ""
    row.update(
        dnm_strict=ev.get("strict") or "",
        DnmStrictBadges=strict_badges(ev.get("strict")),
        dnm_flags=",".join(sorted(flags, key=flag_rank)),
        DnmFlagBadges=flag_badges(flags),
        _dnm_n_flags=len([f for f in flags if f not in STRICT_IGNORED]),
        dnm_child_reads=f"{alt}/{alt + ref}"
        if alt is not None and ref is not None
        else None,
        dnm_child_vaf=_float(ev.get("child_vaf")),
        dnm_father_alt=_int(ev.get("father_alt")),
        dnm_mother_alt=_int(ev.get("mother_alt")),
        dnm_sibs=sibs,
        dnm_origin=ev.get("origin") or "",
        DnmOriginBadges=origin_badges(
            ev.get("origin"), ev.get("origin_status"), ev.get("origin_evidence")
        ),
        dnm_recurrence=_int(ev.get("recurrence_families")),
        dnm_cluster_n=_int(ev.get("cluster_n")),
        dnm_nearby_calls=_int(ev.get("nearby_calls")),
        dnm_binom_p=_float(ev.get("child_binom_p")),
        dnm_strand_p=_float(ev.get("strand_p")),
        _dnm=dict(ev),
    )


def add_evidence(rows: List[Dict[str, Any]], data: DenovoData) -> int:
    """Set every row's evidence; returns how many the step has no row for."""
    by_key = {candidate_key(r): r for r in data.rows()}
    missing = 0
    for row in rows:
        ev = by_key.get((row.get("Variant") or "", str(row.get("sample") or "")))
        if ev is None:
            missing += 1
        set_evidence(row, ev)
    return missing


def row_flags(row: Dict[str, Any]) -> Set[str]:
    return set(split_flags(row.get("dnm_flags")))


@dataclass
class DenovoFilter:
    """The de novo sub-tab's own filter: ``{"value": ...}`` dicts, as menus take.

    ``mode`` is ``strict`` (the default) or ``all``; ``hide`` drops a row with
    any of its flags, ``only`` keeps a row with at least one of its flags.
    """

    mode: Dict[str, str] = field(default_factory=lambda: {"value": "strict"})
    hide: Dict[str, List[str]] = field(default_factory=lambda: {"value": []})
    only: Dict[str, List[str]] = field(default_factory=lambda: {"value": []})

    @property
    def is_default(self) -> bool:
        return (
            self.mode["value"] == "strict"
            and not self.hide["value"]
            and not self.only["value"]
        )


def passes(row: Dict[str, Any], f: DenovoFilter) -> bool:
    """Whether a row passes the de novo filter. A row not assessed is strict."""
    if f.mode["value"] == "strict" and row.get("dnm_strict") == "no":
        return False
    flags = row_flags(row)
    if f.hide["value"] and flags & set(f.hide["value"]):
        return False
    if f.only["value"] and not flags & set(f.only["value"]):
        return False
    return True


def apply_denovo_filter(
    rows: List[Dict[str, Any]], f: DenovoFilter
) -> List[Dict[str, Any]]:
    return [r for r in rows if passes(r, f)]


def is_curated_denovo(row: Dict[str, Any]) -> bool:
    """Validated present without an inherited mode, or carrying a diagnostic.

    What a filter must not hide unnoticed: a curator's de novo call, or one
    behind a diagnosis. An inheritance of ``paternal``, ``maternal``,
    ``either`` or ``homozygous`` says it is not de novo.
    """
    inherited = (row.get("ValidationInheritance") or "") in INHERITED
    if inherited:
        return False
    if row.get("Validation") in ("present", IN_PHASE_MNV):
        return True
    return bool(row.get("Diagnostic"))


def hidden_curated(
    rows: List[Dict[str, Any]], shown: List[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Curated de novo rows (:func:`is_curated_denovo`) a filter left out."""
    ids = {id(r) for r in shown}
    return [r for r in rows if id(r) not in ids and is_curated_denovo(r)]


# ---------------------------------------------------------------------------
# Curation against the candidates
# ---------------------------------------------------------------------------


def _read_records(
    path: Path, families: Optional[Set[str]]
) -> Dict[Tuple[str, str, str], List[Dict[str, str]]]:
    """``(FID, Variant, Sample)`` -> non-ignored rows of a curation file."""
    found: Dict[Tuple[str, str, str], List[Dict[str, str]]] = {}
    if not path.is_file():
        return found
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            fid = row.get("FID") or ""
            if families is not None and fid not in families:
                continue
            if (row.get("Ignore") or "0") == "1":
                continue
            key = (fid, row.get("Variant") or "", row.get("Sample") or "")
            if key[1] and key[2]:
                found.setdefault(key, []).append(row)
    return found


@dataclass
class Curation:
    """The SNV validations (consensus) and diagnostics of some families."""

    validations: Dict[Tuple[str, str, str], Tuple[str, str]]
    diagnostics: Dict[Tuple[str, str, str], str]

    def status(self, key: Tuple[str, str, str]) -> Tuple[str, str]:
        return self.validations.get(key, ("", ""))

    def is_denovo(self, key: Tuple[str, str, str]) -> bool:
        status, inheritance = self.status(key)
        if inheritance in INHERITED:
            return False
        return status in ("present", IN_PHASE_MNV) or key in self.diagnostics


def load_curation(data_dir: Path, families: Optional[Iterable[str]] = None) -> Curation:
    """Blocking: reads ``validations/snvs.tsv`` and ``diagnostics/snvs.tsv``."""
    wanted = set(families) if families is not None else None
    raw = _read_records(data_dir / "validations" / "snvs.tsv", wanted)
    validations = {
        key: consensus(
            (r.get("Validation", ""), r.get("Inheritance", "")) for r in rows
        )
        for key, rows in raw.items()
    }
    diags = _read_records(data_dir / "diagnostics" / "snvs.tsv", wanted)
    diagnostics = {
        key: ", ".join(dict.fromkeys(r.get("Diagnostic", "") for r in rows))
        for key, rows in diags.items()
    }
    return Curation(validations, diagnostics)


def curation_check(
    candidates: List[Dict[str, Any]],
    curation: Curation,
    families: Iterable[str],
) -> List[Dict[str, Any]]:
    """Curated de novo calls the strict set loses, or the candidates lack.

    - ``hidden``: a candidate that is not strict, though validated present
      without an inherited mode or carrying a diagnostic;
    - ``missing``: a validation saying ``de novo`` (present) in a family the
      step assessed, with no candidate at all: pywombat's rule or a newer
      table dropped it.
    """
    families = set(families)
    by_key = {(str(r.get("family_id")), *candidate_key(r)): r for r in candidates}
    out: List[Dict[str, Any]] = []
    keys = set(curation.validations) | set(curation.diagnostics)
    for key in sorted(keys):
        fid, variant, sample = key
        if fid not in families:
            continue
        status, inheritance = curation.status(key)
        cand = by_key.get(key)
        if cand is not None:
            if cand.get("strict") != "no" or not curation.is_denovo(key):
                continue
            problem = "hidden"
        elif status in ("present", IN_PHASE_MNV) and inheritance == "de novo":
            problem = "missing"
        else:
            continue
        out.append(
            {
                "id": "|".join(key),
                "FID": fid,
                "sample": sample,
                "Variant": variant,
                "problem": problem,
                "validation": status,
                "inheritance": inheritance,
                "diagnostic": curation.diagnostics.get(key, ""),
                "flags": (cand or {}).get("flags") or "",
                "FlagBadges": flag_badges(split_flags((cand or {}).get("flags"))),
                "genes": (cand or {}).get("genes") or "",
            }
        )
    return out


# ---------------------------------------------------------------------------
# Cohort summaries
# ---------------------------------------------------------------------------


def _is_transition(ref: str, alt: str) -> bool:
    return {ref, alt} in ({"A", "G"}, {"C", "T"})


def _share(n: int, d: int) -> Optional[float]:
    return round(n / d, 4) if d else None


def child_rows(
    data: DenovoData, keys: Optional[Set[Tuple[str, str]]] = None
) -> List[Dict[str, Any]]:
    """One row per child of the QC table, recounted from the candidates.

    Strict counts are autosomal, as the pipeline's QC; SNV and indel rates are
    per callable unmasked Gb. ``keys`` ``(family_id, child)`` restricts them.
    """
    by_child: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}
    for r in data.rows():
        by_child.setdefault((str(r.get("family_id")), str(r.get("child"))), []).append(
            r
        )
    out = []
    for q in data.qc_rows():
        key = (str(q.get("family_id")), str(q.get("child")))
        if keys is not None and key not in keys:
            continue
        recs = by_child.get(key, [])
        auto = [r for r in recs if r.get("chrom") in AUTOSOMES]
        strict = [r for r in auto if r.get("strict") == "yes"]
        snv = [r for r in strict if r.get("var_type") == "SNV"]
        indel = [r for r in strict if r.get("var_type") == "INDEL"]
        ti = sum(1 for r in snv if _is_transition(str(r.get("ref")), str(r.get("alt"))))
        phased = [r for r in strict if r.get("origin") in ("paternal", "maternal")]
        paternal = sum(1 for r in phased if r.get("origin") == "paternal")
        gb = (_float(q.get("callable_unmasked_bp")) or 0) / 1e9
        flag_sets = [set(split_flags(r.get("flags"))) for r in recs]

        def flagged(flag: str, sets: List[Set[str]] = flag_sets) -> int:
            return sum(1 for f in sets if flag in f)

        info = qc_info(q.get("qc"))
        out.append(
            {
                "id": f"{key[0]}|{key[1]}",
                "FID": key[0],
                "child": key[1],
                "father": q.get("father"),
                "mother": q.get("mother"),
                "sex": q.get("child_sex"),
                "n_candidates": len(recs),
                "n_strict": len(strict),
                "n_strict_snv": len(snv),
                "n_strict_indel": len(indel),
                "rate_per_gb": round(len(strict) / gb, 2) if gb else None,
                "snv_per_gb": round(len(snv) / gb, 2) if gb else None,
                "indel_per_gb": round(len(indel) / gb, 2) if gb else None,
                "indel_share": _share(len(indel), len(strict)),
                "ti_fraction": _share(ti, len(snv)),
                "n_phased": len(phased),
                "phased_share": _share(len(phased), len(strict)),
                "paternal_share": _share(paternal, len(phased)),
                "n_parent_alt": flagged("PARENT_ALT_POSSIBLE")
                + flagged("PARENT_ALT_LIKELY"),
                "n_child_mosaic": flagged("CHILD_MOSAIC_POSSIBLE"),
                "n_sib_shared": flagged("SIB_SHARED"),
                "n_sib_shared_strict": sum(
                    1 for r in strict if "SIB_SHARED" in split_flags(r.get("flags"))
                ),
                "callable_gb": round(gb, 3) if gb else None,
                "callable_status": q.get("callable_status"),
                "qc": q.get("qc") or "",
                "QcBadges": [info] if q.get("qc") else [],
                "qc_low": _float(q.get("qc_low")),
                "qc_high": _float(q.get("qc_high")),
                "qc_reference": q.get("qc_reference"),
                "sampleqc_flags": q.get("sampleqc_flags") or "",
            }
        )
    return out


def _median(values: Iterable[Optional[float]]) -> Optional[float]:
    vals = [v for v in values if v is not None]
    return statistics.median(vals) if vals else None


def cohort_summary(children: List[Dict[str, Any]]) -> Dict[str, Optional[float]]:
    """Medians over children and pooled shares, for the literature table."""
    strict = sum(c["n_strict"] for c in children)
    snv = sum(c["n_strict_snv"] for c in children)
    indel = sum(c["n_strict_indel"] for c in children)
    phased = sum(c["n_phased"] for c in children)
    paternal = sum(round((c["paternal_share"] or 0) * c["n_phased"]) for c in children)
    return {
        "children": len(children),
        "strict_median": _median(c["n_strict"] for c in children),
        "snv_median": _median(c["n_strict_snv"] for c in children),
        "indel_median": _median(c["n_strict_indel"] for c in children),
        "rate_median": _median(c["rate_per_gb"] for c in children),
        "snv_rate_median": _median(c["snv_per_gb"] for c in children),
        "indel_rate_median": _median(c["indel_per_gb"] for c in children),
        "indel_share": _share(indel, strict),
        "ti_fraction": _median(c["ti_fraction"] for c in children),
        "phased_share": _share(phased, strict),
        "paternal_share": _share(paternal, phased),
        "callable_gb_median": _median(c["callable_gb"] for c in children),
        "implied_snv_mu": (
            snv_rate / 2e9
            if (snv_rate := _median(c["snv_per_gb"] for c in children))
            else None
        ),
        "n_strict": strict,
        "n_snv": snv,
        "n_indel": indel,
        "parent_mosaic": _share(
            sum(c["n_parent_alt"] for c in children),
            sum(c["n_candidates"] for c in children),
        ),
        "child_mosaic": _share(
            sum(c["n_child_mosaic"] for c in children),
            sum(c["n_candidates"] for c in children),
        ),
        "sib_shared": _share(sum(c["n_sib_shared_strict"] for c in children), strict),
    }


def flag_summary(
    candidates: List[Dict[str, Any]], curation: Optional[Curation] = None
) -> List[Dict[str, Any]]:
    """Per flag: candidates it marks, alone, and what curators said of them."""
    total = len(candidates)
    out: Dict[str, Dict[str, Any]] = {}
    for r in candidates:
        flags = split_flags(r.get("flags"))
        key = (str(r.get("family_id")), *candidate_key(r))
        status, inheritance = curation.status(key) if curation else ("", "")
        present = status in ("present", IN_PHASE_MNV) and inheritance not in INHERITED
        for flag in flags:
            entry = out.setdefault(
                flag,
                {
                    "flag": flag,
                    "n": 0,
                    "alone": 0,
                    "present": 0,
                    "absent": 0,
                    "diagnosed": 0,
                },
            )
            entry["n"] += 1
            if len(flags) == 1:
                entry["alone"] += 1
            if present:
                entry["present"] += 1
            elif status == "absent":
                entry["absent"] += 1
            if curation and key in curation.diagnostics:
                entry["diagnosed"] += 1
    rows = sorted(out.values(), key=lambda e: flag_rank(e["flag"]))
    for e in rows:
        info = flag_info(e["flag"])
        e.update(
            id=e["flag"],
            label=info["label"],
            share=_share(e["n"], total),
            FlagBadges=flag_badges([e["flag"]]),
            category=CATEGORIES.get(info["category"], {}).get(
                "label", info["category"]
            ),
        )
    return rows


def strict_vs_curation(
    candidates: List[Dict[str, Any]], curation: Curation
) -> Dict[str, Dict[str, int]]:
    """Strict / flagged against validated present (de novo) / absent."""
    table = {s: {"present": 0, "absent": 0, "other": 0} for s in ("yes", "no")}
    for r in candidates:
        key = (str(r.get("family_id")), *candidate_key(r))
        if key not in curation.validations:
            continue
        status, inheritance = curation.status(key)
        col = (
            "present"
            if status in ("present", IN_PHASE_MNV) and inheritance not in INHERITED
            else "absent"
            if status == "absent"
            else "other"
        )
        table.setdefault(
            r.get("strict") or "no", {"present": 0, "absent": 0, "other": 0}
        )[col] += 1
    return table


def candidates_export(rows: List[Dict[str, Any]], columns: List[str]) -> pl.DataFrame:
    """The candidates as a frame with ``columns`` (the ones the table has)."""
    if not rows:
        return pl.DataFrame({c: [] for c in columns})
    df = pl.DataFrame(rows, infer_schema_length=None)
    return df.select([c for c in columns if c in df.columns])
