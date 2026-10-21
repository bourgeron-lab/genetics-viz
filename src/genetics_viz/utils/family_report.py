"""The family report: what the family page shows, gathered for a clinician.

The report goes to the referring clinician of a family's centre, who decides
whether to recontact the family. It is not addressed to patients. It recaps:
- the pedigree;
- ancestry;
- Sample QC issues;
- runs of homozygosity;
- a few polygenic scores;
- aneuploidies;
- UPDs;
- the curators' pathogenic and uncertain diagnostics;
- short methods.

Only curators produce it, so it uses the curator wording throughout.

:func:`build_family_report` is blocking and pure. It reads no
``app.storage``: the caller resolves the user on the event loop and runs it
in ``asyncio.to_thread``. Each section is a :class:`Section`, either
``available`` with rows, or not, with the reason. A missing file is a
degraded state, as in the tabs, never an error. :mod:`report_html` renders
the result.

Wording, the traits shown and the methods live in ``config/report.yaml``.
"""

import logging
import math
from dataclasses import dataclass, field
from datetime import datetime
from importlib import metadata
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import polars as pl
import yaml

from genetics_viz.components.diagnostic_loader import load_family_diagnostics
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)

from . import roh as roh_utils
from . import sample_qc as qc_utils
from . import upd as upd_utils
from .ancestry import (
    AncestryFiles,
    ancestry_export,
    find_ancestry_files,
    get_region_color,
    load_reference_pcs,
    pgs_export,
)
from .ancestry_plots import zoom_window
from .cohort_state import read_cohort_state
from .cytobands import CHROM_ORDER, CYTOBANDS, get_cytoband_range, norm_chrom
from .diagnostic_labels import describe_str, describe_upd
from .export import person
from .pedigree_labels import PHENO_COLORS, normalize_pheno, normalize_sex
from .pipeline_params import load_params
from .report_summaries import FamilySummary, latest_family_summary
from .report_svg import Marker, ordinal, pca_svg, percentile, pgs_distribution_svg
from .sample_qc_acks import ack_summary, load_acks, split_acknowledged
from .sharding import get_family_path
from .str_eh import find_eh_labels, get_family_eh_dir
from .validation_files import validation_path
from .wisecondorx import parse_wisecondorx_bed

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "report.yaml"


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


CONFIG: Dict[str, Any] = _load_config()

_SEX_TEXT = {"1": "Male", "2": "Female", "-9": "Unknown"}
_PHENO_TEXT = {"2": "Affected", "1": "Unaffected", "-9": "Unknown"}

#: Acrocentric chromosomes: their p arm holds no assessable sequence.
_ACROCENTRIC = {"13", "14", "15", "21", "22"}


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


@dataclass
class Plot:
    caption: str
    svg: str
    note: str = ""


@dataclass
class Section:
    """One block of the report."""

    key: str
    title: str
    available: bool = True
    #: Why the section is empty or missing, or what it found in one sentence.
    note: str = ""
    #: ``(key, header)`` of the table, in order.
    columns: List[Tuple[str, str]] = field(default_factory=list)
    rows: List[Dict[str, Any]] = field(default_factory=list)
    plots: List[Plot] = field(default_factory=list)
    #: Short paragraphs under the table.
    footnotes: List[str] = field(default_factory=list)
    #: Something the clinician should look at (drawn as an alert).
    attention: bool = False
    #: The line under "Results at a glance"; the note when empty.
    glance: str = ""

    @classmethod
    def missing(cls, key: str, title: str, reason: str) -> "Section":
        return cls(key=key, title=title, available=False, note=reason)


@dataclass
class FamilyReport:
    family_id: str
    cohort_name: str
    generated_at: str
    generated_by: str
    app_version: str
    pipeline: str
    summary: Optional[FamilySummary]
    sections: List[Section]
    methods: List[Tuple[str, str]]
    title: str = CONFIG.get("title", "")
    subtitle: str = CONFIG.get("subtitle", "")
    audience: str = CONFIG.get("audience", "")
    disclaimer: str = CONFIG.get("disclaimer", "")
    methods_source: str = CONFIG.get("methods_source", "")

    @property
    def document_title(self) -> str:
        """The browser's default PDF file name."""
        parts = [self.cohort_name, self.family_id, "report", self.generated_at[:10]]
        return "_".join(p.replace(" ", "-") for p in parts if p)

    def section(self, key: str) -> Optional[Section]:
        return next((s for s in self.sections if s.key == key), None)


@dataclass
class _Member:
    iid: str
    father: str
    mother: str
    sex: str
    phenotype: str

    @property
    def color(self) -> str:
        return PHENO_COLORS.get(self.phenotype, PHENO_COLORS["-9"])


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _float(value: Any) -> Optional[float]:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _fmt(value: Any, digits: int = 3) -> str:
    f = _float(value)
    return "—" if f is None else f"{f:.{digits}f}"


def _truthy(value: Any) -> bool:
    """``is_outlier`` is written as ``True``/``False`` text by the pipeline."""
    if isinstance(value, str):
        return value.strip().lower() in ("true", "1", "yes")
    return bool(value)


def _by_value(rows: List[Dict[str, Any]], key: str) -> str:
    """``"EUR: S1, S2 · CSA: S3"``, or ``"EUR (3 individuals)"`` for one value."""
    groups: Dict[str, List[str]] = {}
    for r in rows:
        groups.setdefault(str(r.get(key) or "—"), []).append(str(r.get("iid") or ""))
    if len(groups) == 1:
        ((value, ids),) = groups.items()
        return f"{value} ({len(ids)} individual{'s' if len(ids) > 1 else ''})"
    return " · ".join(f"{value}: {', '.join(ids)}" for value, ids in groups.items())


def _app_version() -> str:
    try:
        return metadata.version("genetics-viz")
    except metadata.PackageNotFoundError:
        return "unknown"


def _pipeline_text(cohort: Any) -> str:
    if cohort is None:
        return "not recorded (family opened outside its cohort)"
    try:
        state = read_cohort_state(cohort.path)
    except Exception:
        logger.warning(
            "Failed to read the run record of %s", cohort.name, exc_info=True
        )
        state = None
    run = state.reference_run if state is not None and not state.unsupported else None
    if run is None:
        return "not recorded"
    text = f"ghfc-ngs {run.pipeline_version or '?'}"
    if run.short_commit:
        text += f" (commit {run.short_commit})"
    if run.finished_at is not None:
        text += f", run of {run.finished_at.date().isoformat()}"
    return text


class _Defaults(dict):
    def __missing__(self, key: str) -> str:
        return "n/a"


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------


def _pedigree(members: List[_Member], qc: Optional[qc_utils.QcData]) -> Section:
    karyotype = {r.get("sample"): r for r in (qc.sample_rows() if qc else [])}
    rows = []
    for m in members:
        k = karyotype.get(m.iid) or {}
        rows.append(
            {
                "iid": m.iid,
                "father": m.father or "—",
                "mother": m.mother or "—",
                "sex": _SEX_TEXT[m.sex],
                "phenotype": _PHENO_TEXT[m.phenotype],
                "_phenotype_color": m.color,
                "karyotype": k.get("karyotype") or "—",
            }
        )
    return Section(
        key="pedigree",
        title="Family",
        columns=[
            ("iid", "Individual"),
            ("father", "Father"),
            ("mother", "Mother"),
            ("sex", "Sex (declared)"),
            ("phenotype", "Status"),
            ("karyotype", "Sex chromosomes (observed)"),
        ],
        rows=rows,
    )


def _people(members: List[_Member], family_id: str) -> List[Dict[str, Any]]:
    return [person(family_id, m.iid, m.sex, m.phenotype) for m in members]


def _ancestry(
    files: Optional[AncestryFiles], members: List[_Member], family_id: str
) -> Section:
    title = "Genetic ancestry"
    if files is None or files.pcs_path is None:
        return Section.missing(
            "ancestry", title, "No ancestry results for this family."
        )
    try:
        df = ancestry_export(files, _people(members, family_id))
    except Exception:
        logger.warning("Failed to read the ancestry of %s", family_id, exc_info=True)
        return Section.missing(
            "ancestry", title, "The ancestry results could not be read."
        )
    admix_min = float(CONFIG.get("admixture_min", 0.05))
    admix_cols = [c for c in df.columns if c.startswith("admix_")]
    rows = []
    for r in df.to_dicts():
        parts = sorted(
            ((c[len("admix_") :], _float(r.get(c))) for c in admix_cols),
            key=lambda p: -(p[1] or 0),
        )
        admixture = " · ".join(
            f"{name} {v * 100:.0f}%"
            for name, v in parts
            if v is not None and v >= admix_min
        )
        rows.append(
            {
                "iid": r["IID"],
                "region": r.get("region") or "—",
                "_region_color": get_region_color(r.get("region")),
                "population": r.get("population") or "—",
                "confidence": _fmt(r.get("region_confidence"), 2),
                "admixture": admixture or "—",
                "outlier": "Yes" if _truthy(r.get("is_outlier")) else "",
            }
        )
    section = Section(
        key="ancestry",
        title=title,
        columns=[
            ("iid", "Individual"),
            ("region", "Region"),
            ("population", "Closest population"),
            ("confidence", "Confidence"),
            ("admixture", "Admixture"),
            ("outlier", "Outlier"),
        ],
        rows=rows,
        note=f"Reference bundle {files.label}.",
        glance=_by_value(rows, "region"),
    )
    outliers = [r["iid"] for r in rows if r["outlier"]]
    if outliers:
        section.footnotes.append(
            f"{', '.join(outliers)}: far from every reference genome; the region is"
            " the nearest one, not a close match."
        )
    section.plots = _pca_plots(df, files, rows)
    return section


def _pca_plots(df: pl.DataFrame, files: AncestryFiles, rows: List[Dict]) -> List[Plot]:
    if not {"PC1", "PC2"} <= set(df.columns):
        return []
    points = df.select(["IID", "PC1", "PC2"]).to_dicts()
    colors = {r["iid"]: r["_region_color"] for r in rows}
    markers = [
        (p["IID"], x, y, colors.get(p["IID"], "#111827"))
        for p in points
        if (x := _float(p["PC1"])) is not None and (y := _float(p["PC2"])) is not None
    ]
    if not markers:
        return []
    reference = load_reference_pcs(files.bundle_version)
    ref_points = None
    regions = sorted({r["region"] for r in rows if r["region"] != "—"})
    window = None
    note = ""
    if reference is not None:
        ref_points = [
            (x, y, str(region))
            for x, y, region in zip(
                reference["PC1"].to_list(),
                reference["PC2"].to_list(),
                reference["region"].to_list(),
            )
            if x is not None and y is not None
        ]
        window = zoom_window(points, reference, regions, "PC1", "PC2")
    else:
        note = "Reference panel unavailable: the family is drawn alone."
    plots = []
    if ref_points:
        plots.append(Plot("Whole reference panel", pca_svg(markers, ref_points, None)))
    plots.append(
        Plot(
            "Zoomed on the family's region" if ref_points else "Family",
            pca_svg(markers, ref_points, window),
            note,
        )
    )
    return plots


def _sample_qc(
    data_dir: Path,
    family_id: str,
    family_qc: Optional[qc_utils.QcData],
    cohort: Any,
) -> Section:
    title = "Sample quality control"
    if family_qc is None:
        return Section.missing(
            "sample_qc", title, "No sample QC results for this family."
        )
    cohort_qc = None
    if cohort is not None:
        try:
            cohort_qc = qc_utils.load_cohort_qc(cohort)
        except Exception:
            logger.warning(
                "Failed to read the cohort QC of %s", cohort.name, exc_info=True
            )
    flags = qc_utils.family_flags_with_cohort(family_qc, cohort_qc, family_id)
    acks = load_acks(data_dir)
    active, acked = split_acknowledged(flags, acks)
    rows = []
    for f in active:
        sev = qc_utils.severity_info(f.severity)
        rows.append(
            {
                "severity": sev["label"],
                "_severity_color": sev["color"],
                "samples": " – ".join(f.samples),
                "issue": qc_utils.flag_title(f, True),
                "_issue_note": f.info.explain,
                "detail": f.evidence,
                "review": "Not reviewed",
            }
        )
    for f in acked:
        rows.append(
            {
                "severity": "Reviewed",
                "_severity_color": "#9e9e9e",
                "samples": " – ".join(f.samples),
                "issue": qc_utils.flag_title(f, True),
                "detail": f.evidence,
                "review": ack_summary(f, acks[f], True),
                "_muted": True,
            }
        )
    counts = qc_utils.count_by_severity(active)
    summary = qc_utils.severity_summary(counts)
    if not rows:
        note = CONFIG["none"]["sample_qc"]
    elif summary:
        note = f"Open issues: {summary}."
    else:
        note = "Every issue raised has been reviewed by a curator."
    if acked:
        note += f" {len(acked)} reviewed and explained or confirmed."
    if cohort is None:
        note += " Relatedness with other families was not checked (no cohort)."
    section = Section(
        key="sample_qc",
        title=title,
        note=note,
        columns=[
            ("severity", "Severity"),
            ("samples", "Individual(s)"),
            ("issue", "Issue"),
            ("detail", "Evidence"),
            ("review", "Curator review"),
        ],
        rows=rows,
        attention=bool(counts.get("high")),
    )
    metrics = []
    for r in family_qc.sample_rows():
        metrics.append(
            {
                "iid": r.get("sample"),
                "depth": _fmt(r.get("depth_mean"), 1),
                "contamination": _fmt(r.get("contamination_charr"), 3),
                "mendelian": _fmt(r.get("trio_mendelian_error_rate"), 4),
            }
        )
    if metrics:
        section.footnotes.append(
            "Per individual — mean depth (somalier sites) / contamination (CHARR) /"
            " trio Mendelian error rate: "
            + "; ".join(
                f"{m['iid']}: {m['depth']}× / {m['contamination']} / {m['mendelian']}"
                for m in metrics
            )
            + "."
        )
    return section


def _roh(
    data_dir: Path, members: List[_Member], family_id: str
) -> Tuple[Section, Optional[pl.DataFrame]]:
    title = "Runs of homozygosity"
    files = roh_utils.find_roh_files(data_dir, family_id)
    froh = None
    if files is not None:
        try:
            froh = roh_utils.read_roh_table(files.path("froh"))
        except Exception:
            logger.warning("Failed to read the ROH of %s", family_id, exc_info=True)
    if froh is None or len(froh) == 0:
        return Section.missing("roh", title, "No ROH results for this family."), None
    people = {p["IID"]: p for p in _people(members, family_id)}
    paired = [(r, people[r["IID"]]) for r in froh.to_dicts() if r.get("IID") in people]
    export_rows, _ = roh_utils.froh_export(paired, show_relationship=True)
    rows = []
    for r in export_rows:
        grade = roh_utils.grade_for(r.get("froh_category"))
        pct = _float(r.get("pct_froh_1.5mb_vs_ref_group"))
        rows.append(
            {
                "iid": r["IID"],
                "froh": _fmt(r.get("froh_1.5mb"), 3),
                "pct": (
                    f"{roh_utils.ordinal(pct)} (vs {r.get('ref_group') or '?'},"
                    f" n={r.get('ref_n') or '?'})"
                    if pct is not None
                    else "—"
                ),
                "pattern": r.get("pattern") or "—",
                "grade": grade.grade,
                "_grade_color": grade.color,
                "_grade_note": (
                    grade.relationship
                    if r.get("froh_category") not in (None, "none")
                    else ""
                ),
                "kinship": (
                    f"{_fmt(r.get('parents_kinship'), 3)} (expected F"
                    f" {_fmt(r.get('expected_f'), 3)})"
                    if _float(r.get("parents_kinship")) is not None
                    else "—"
                ),
                "flags": r.get("flags") or "",
            }
        )
    flagged = [r for r in paired if roh_utils.is_flagged(r[0])]
    note = (
        f"Unusual homozygosity in {len(flagged)} of {len(rows)} individual(s)."
        if flagged
        else "Homozygosity within the range of the reference population for everyone."
    )
    return (
        Section(
            key="roh",
            title=title,
            note=note,
            columns=[
                ("iid", "Individual"),
                ("froh", "F_ROH"),
                ("pct", "Percentile"),
                ("pattern", "Pattern"),
                ("grade", "Grade"),
                ("kinship", "Parents' kinship"),
                ("flags", "Flags"),
            ],
            rows=rows,
            footnotes=[roh_utils.GRADE_HELP] if roh_utils.GRADE_HELP else [],
        ),
        froh,
    )


def _pgs(
    files: Optional[AncestryFiles], members: List[_Member], family_id: str
) -> Section:
    title = "Polygenic scores"
    traits = CONFIG.get("pgs_traits") or []
    if files is None or files.pgs_zscore_path is None:
        return Section.missing("pgs", title, "No polygenic scores for this family.")
    try:
        df = pgs_export(files, _people(members, family_id))
    except Exception:
        logger.warning("Failed to read the PGS of %s", family_id, exc_info=True)
        return Section.missing("pgs", title, "The polygenic scores could not be read.")
    by_iid = {r["IID"]: r for r in df.to_dicts()}
    colors = {m.iid: m.color for m in members}
    rows = [{"iid": m.iid} for m in members if m.iid in by_iid]
    columns = [("iid", "Individual")]
    plots = []
    missing = []
    extremes: List[str] = []
    for t in traits:
        tid, label = t["id"], t.get("label", t["id"])
        if tid not in df.columns:
            missing.append(label)
            continue
        columns.append((tid, label))
        markers = []
        for row in rows:
            z = _float(by_iid[row["iid"]].get(tid))
            if z is None:
                row[tid] = "—"
                continue
            row[tid] = f"z {z:+.2f} · {ordinal(percentile(z))}"
            if abs(z) >= 2:
                extremes.append(f"{row['iid']} {label} {ordinal(percentile(z))}")
            markers.append(
                Marker(f"{row['iid']} {ordinal(percentile(z))}", z, colors[row["iid"]])
            )
        plots.append(Plot(label, pgs_distribution_svg(markers), t.get("source", "")))
    if len(columns) == 1:
        return Section.missing(
            "pgs", title, "None of the selected scores is in this family's results."
        )
    footnotes = [CONFIG.get("pgs_note", "")]
    if missing:
        footnotes.append(f"Not computed for this family: {', '.join(missing)}.")
    return Section(
        key="pgs",
        title=title,
        note="Markers: red affected, green unaffected, grey unknown.",
        glance="; ".join(extremes)
        or "No score beyond two standard deviations of the reference.",
        columns=columns,
        rows=rows,
        plots=plots,
        footnotes=[f for f in footnotes if f],
    )


def _assessable_span(chrom: str) -> Optional[Tuple[float, float]]:
    """``(start, end)`` in Mb of what read depth can assess: the q arm of an
    acrocentric chromosome, the whole chromosome otherwise."""
    bands = CYTOBANDS.get(chrom) or []
    if not bands:
        return None
    end = max(b["end"] for b in bands)
    start = 0.0
    if chrom in _ACROCENTRIC:
        start = max((b["end"] for b in bands if b["stain"] == "acen"), default=0.0)
    return start, end


def whole_chromosome_events(
    calls: Optional[pl.DataFrame], min_fraction: float
) -> Dict[str, List[str]]:
    """``sample -> ["+21", "−18"]``: gains/losses covering most of a chromosome.

    Same-direction segments (gain for a positive log2 ratio) are summed per
    sample and chromosome, against the chromosome's assessable length.
    """
    out: Dict[str, List[str]] = {}
    if calls is None or len(calls) == 0:
        return out
    needed = {"chr", "start", "end", "ratio", "sample"}
    if not needed <= set(calls.columns):
        return out
    covered: Dict[Tuple[str, str, str], float] = {}
    for r in calls.to_dicts():
        ratio, start, end = _float(r["ratio"]), _float(r["start"]), _float(r["end"])
        if ratio is None or start is None or end is None or ratio == 0:
            continue
        chrom = norm_chrom(r["chr"])
        span = _assessable_span(chrom)
        if span is None:
            continue
        overlap = min(end / 1e6, span[1]) - max(start / 1e6, span[0])
        if overlap > 0:
            key = (str(r["sample"]), chrom, "gain" if ratio > 0 else "loss")
            covered[key] = covered.get(key, 0.0) + overlap
    order = {c: i for i, c in enumerate(CHROM_ORDER)}
    for (sample, chrom, kind), mb in sorted(
        covered.items(), key=lambda kv: (kv[0][0], order.get(kv[0][1], 99))
    ):
        span = _assessable_span(chrom)
        length = span[1] - span[0] if span else 0.0
        if length > 0 and mb / length >= min_fraction:
            out.setdefault(sample, []).append(
                f"{'+' if kind == 'gain' else '−'}{chrom} ({mb / length * 100:.0f}%)"
            )
    return out


def _aneuploidy(
    data_dir: Path,
    members: List[_Member],
    family_id: str,
    family_qc: Optional[qc_utils.QcData],
) -> Section:
    title = "Aneuploidy"
    wcx = get_family_path(data_dir, family_id) / "svs" / "wisecondorx"
    wcx_file = wcx / f"{family_id}_aberrations.annotated.bed"
    calls = None
    if wcx_file.exists():
        try:
            calls = parse_wisecondorx_bed(wcx_file)
        except Exception:
            logger.warning("Failed to read %s", wcx_file, exc_info=True)
    if family_qc is None and not wcx_file.exists():
        return Section.missing(
            "aneuploidy", title, "No copy-number or sample QC results for this family."
        )
    min_fraction = float(CONFIG.get("aneuploidy_min_fraction", 0.75))
    events = whole_chromosome_events(calls, min_fraction)
    qc_rows = {
        r.get("sample"): r for r in (family_qc.sample_rows() if family_qc else [])
    }
    rows = []
    #: One "S5: XXY, +21" per individual with something to report.
    glance = []
    for m in members:
        q = qc_rows.get(m.iid)
        found_here = []
        if q is None:
            sex_chroms = "not assessed"
        else:
            karyotype = q.get("karyotype") or "unknown"
            expected = {"1": "XY", "2": "XX"}.get(m.sex)
            sex_chroms = (
                f"{karyotype} (X {_fmt(q.get('x_copies'), 2)},"
                f" Y {_fmt(q.get('y_copies'), 2)})"
            )
            if not qc_utils.is_typical_karyotype(karyotype):
                found_here.append(karyotype)
            elif expected and karyotype != expected:
                sex_chroms += f" — declared {_SEX_TEXT[m.sex].lower()}"
                found_here.append(f"{karyotype}, declared {_SEX_TEXT[m.sex].lower()}")
        if not wcx_file.exists():
            autosomes = "not assessed"
        else:
            autosomes = ", ".join(events.get(m.iid, [])) or "none detected"
            found_here += events.get(m.iid, [])
        if found_here:
            glance.append(f"{m.iid}: {', '.join(found_here)}")
        rows.append(
            {
                "iid": m.iid,
                "sex_chroms": sex_chroms,
                "autosomes": autosomes,
                "_attention": bool(found_here),
            }
        )
    found = bool(glance)
    return Section(
        key="aneuploidy",
        title=title,
        note="" if found else CONFIG["none"]["aneuploidy"],
        columns=[
            ("iid", "Individual"),
            ("sex_chroms", "Sex chromosomes (copies)"),
            ("autosomes", "Whole-chromosome gains / losses"),
        ],
        rows=rows,
        footnotes=[CONFIG.get("aneuploidy_note", "")],
        attention=found,
        glance="; ".join(glance),
    )


def _upd(
    data_dir: Path,
    family_id: str,
    members: List[_Member],
    froh: Optional[pl.DataFrame],
    family_qc: Optional[qc_utils.QcData],
    diagnostics: List[Dict[str, Any]],
) -> Section:
    title = "Uniparental disomy"
    if froh is None and family_qc is None:
        return Section.missing(
            "upd", title, "No ROH or sample QC results for this family."
        )
    member_ids = {m.iid for m in members}
    candidates: Dict[Tuple[str, str], None] = {}
    for r in froh.to_dicts() if froh is not None else []:
        if r.get("upd_candidate_chrom") and r.get("IID") in member_ids:
            candidates[
                (r["IID"], upd_utils.normalise_chrom(r["upd_candidate_chrom"]))
            ] = None
    for f in family_qc.flags if family_qc is not None else []:
        if upd_utils.is_upd_flag(f):
            chrom = upd_utils.chrom_from_evidence(f)
            if chrom:
                candidates[(upd_utils.child_of(f), chrom)] = None
    validations = upd_utils.read_upd_validations(data_dir, family_id)
    for v in validations:
        if v.get("Sample") in member_ids and v.get("Ignore", "0") != "1":
            candidates[
                (v["Sample"], upd_utils.chrom_of_key(v.get("Variant") or ""))
            ] = None
    upd_diag = {
        (s.strip(), e.get("Variant")): e.get("Diagnostic", "")
        for e in diagnostics
        if e.get("_source") == "upd"
        for s in (e.get("Sample") or "").split(",")
    }
    rows = []
    for child, chrom in candidates:
        if not chrom:
            continue
        evidence = upd_utils.gather_evidence(data_dir, family_id, child, chrom)
        key = upd_utils.upd_key(chrom)
        latest = upd_utils.latest_validation(
            [v for v in validations if v.get("Variant") == key], child
        )
        validation = (latest or {}).get("Validation") or "not reviewed"
        rows.append(
            {
                "iid": child,
                "chrom": f"chr{chrom}",
                "evidence": evidence.source,
                "validation": validation,
                "origin": upd_utils.origin_label(
                    (latest or {}).get("Inheritance") or ""
                )
                or "—",
                "diagnostic": upd_diag.get((child, key)) or "—",
                "_muted": validation == "absent",
            }
        )
    member_order = {m.iid: i for i, m in enumerate(members)}
    chrom_order = {c: i for i, c in enumerate(CHROM_ORDER)}
    rows.sort(
        key=lambda r: (
            member_order.get(r["iid"], len(member_order)),
            chrom_order.get(r["chrom"][3:], len(chrom_order)),
        )
    )
    return Section(
        key="upd",
        title=title,
        note="" if rows else CONFIG["none"]["upd"],
        columns=[
            ("iid", "Individual"),
            ("chrom", "Chromosome"),
            ("evidence", "Evidence"),
            ("validation", "Curator validation"),
            ("origin", "Parental origin"),
            ("diagnostic", "Diagnostic"),
        ],
        rows=rows,
        attention=any(r["validation"] in ("present", "not reviewed") for r in rows),
        glance="; ".join(f"{r['iid']} {r['chrom']} ({r['validation']})" for r in rows),
    )


def _genes(gene_str: str, limit: int = 8) -> str:
    symbols = [
        g.split(":")[0].strip()
        for g in (gene_str or "").replace(";", ",").split(",")
        if g.split(":")[0].strip() not in ("", "-", ".")
    ]
    symbols = list(dict.fromkeys(symbols))
    if len(symbols) > limit:
        return ", ".join(symbols[:limit]) + f" (+{len(symbols) - limit})"
    return ", ".join(symbols)


def _sv_text(variant: str) -> str:
    parts = variant.split(":")
    try:
        chrom, span = parts[0], parts[1]
        start, end = (int(x) for x in span.split("-", 1))
    except (IndexError, ValueError):
        return variant
    kind = parts[2].upper() if len(parts) > 2 else "SV"
    kb = (end - start) / 1000
    size = f"{kb / 1000:.2f} Mb" if kb >= 1000 else f"{kb:.0f} kb"
    band = get_cytoband_range(chrom, start, end)
    return f"{kind} {chrom}:{start}-{end} ({', '.join(p for p in (band, size) if p)})"


def _diagnostics(
    data_dir: Path, family_id: str, diagnostics: List[Dict[str, Any]]
) -> Section:
    verdicts: Dict[str, Dict[str, str]] = CONFIG.get("diagnostic_verdicts") or {}
    kept = [e for e in diagnostics if e.get("Diagnostic") in verdicts]
    validation_maps: Dict[str, Dict] = {}
    rows = []
    for e in kept:
        source = e.get("_source", "")
        kind = source.upper()
        variant = e.get("Variant", "")
        if source == "str":
            text, _ = describe_str(data_dir, family_id, e)
            finding = f"{_genes(e.get('Gene', '')) or variant} repeat" + (
                f", {text}" if text else ""
            )
        elif source == "upd":
            finding = f"UPD {describe_upd(data_dir, family_id, e)[0]}"
        elif source == "sv":
            finding = _sv_text(variant)
        else:
            finding = variant
        if kind not in validation_maps:
            try:
                validation_maps[kind] = load_validation_map(
                    validation_path(data_dir, kind), family_id
                )
            except (KeyError, OSError):
                validation_maps[kind] = {}
        validations = []
        for sample in (s.strip() for s in (e.get("Sample") or "").split(",")):
            if not sample:
                continue
            row: Dict[str, Any] = {}
            add_validation_status_to_row(row, validation_maps[kind], variant, sample)
            status = row.get("Validation") or ""
            if status:
                inheritance = row.get("ValidationInheritance") or ""
                validations.append(
                    f"{sample}: {status}" + (f" ({inheritance})" if inheritance else "")
                )
        verdict = e.get("Diagnostic", "")
        rows.append(
            {
                "verdict": verdicts[verdict].get("label", verdict),
                "_verdict_color": verdicts[verdict].get("color", "#757575"),
                "_rank": list(verdicts).index(verdict),
                "kind": kind,
                "finding": finding,
                "genes": _genes(e.get("Gene", "")) if source != "str" else "",
                "impact": e.get("Impact", "") if source in ("snv", "sv") else "",
                "samples": e.get("Sample", ""),
                "validation": "; ".join(validations) or "not validated",
                "curators": e.get("User", ""),
                "comment": e.get("Comment", ""),
            }
        )
    rows.sort(key=lambda r: r["_rank"])
    return Section(
        key="diagnostics",
        title="Findings recorded by the curators",
        note="" if rows else CONFIG["none"]["diagnostics"],
        columns=[
            ("verdict", "Classification"),
            ("kind", "Type"),
            ("finding", "Finding"),
            ("genes", "Gene(s)"),
            ("impact", "Consequence"),
            ("samples", "Individual(s)"),
            ("validation", "Validation (inheritance)"),
            ("comment", "Curator comment"),
        ],
        rows=rows,
        attention=any(r["_rank"] == 0 for r in rows),
        glance=" · ".join(
            f"{n} {label.lower()}"
            for label in dict.fromkeys(r["verdict"] for r in rows)
            if (n := sum(r["verdict"] == label for r in rows))
        ),
    )


def _methods(
    sections: Sequence[Section],
    files: Optional[AncestryFiles],
    family_qc: Optional[qc_utils.QcData],
    data_dir: Path,
    family_id: str,
    cohort: Any,
) -> List[Tuple[str, str]]:
    available = {s.key for s in sections if s.available}
    diagnostics = next((s for s in sections if s.key == "diagnostics"), None)
    kinds = {r["kind"] for r in diagnostics.rows} if diagnostics else set()
    params: Dict[str, Any] = {}
    if cohort is not None and getattr(cohort, "params_file", None):
        params = load_params(cohort.params_file)
    try:
        catalogs = find_eh_labels(get_family_eh_dir(data_dir, family_id), family_id)
    except OSError:
        catalogs = []
    values = _Defaults(
        bundle=files.bundle_version if files else "n/a",
        tag=files.filter_tag if files else "n/a",
        qc_label=family_qc.label if family_qc else "n/a",
        catalog=catalogs[0] if catalogs else "n/a",
        wcx_bins=(
            f"{params['wisecondorx_binsize']} bp bins"
            if params.get("wisecondorx_binsize")
            else "fixed-size bins"
        ),
        wcx_args=params.get("wisecondorx_predict_args", "default settings"),
        aneuploidy_pct=f"{float(CONFIG.get('aneuploidy_min_fraction', 0.75)) * 100:.0f}",
    )
    wanted = ["sequencing"]
    if "aneuploidy" in available:
        wanted.append("cnv")
    for key in ("sample_qc", "ancestry", "roh"):
        if key in available:
            wanted.append(key)
    if "upd" in available:
        wanted.append("upd")
    if "pgs" in available:
        wanted.append("pgs")
    if "STR" in kinds:
        wanted.append("str")
    methods = CONFIG.get("methods") or {}
    out = []
    for key in wanted:
        m = methods.get(key)
        if m:
            out.append((m.get("title", key), str(m.get("text", "")).format_map(values)))
    return out


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def build_family_report(
    data_dir: Path,
    family_id: str,
    family: Any,
    cohort: Any = None,
    user: str = "",
) -> FamilyReport:
    """Every section of the family report. Blocking I/O; never raises on data.

    ``family`` is the :class:`models.Family` (its samples give the members,
    in pedigree order). ``cohort`` is ``None`` on the standalone family page:
    no cross-family QC flags, no run record, no workflow parameters then.
    """
    members = [
        _Member(
            iid=s.sample_id,
            father=s.father_id or "",
            mother=s.mother_id or "",
            sex=normalize_sex(s.sex),
            phenotype=normalize_pheno(s.phenotype),
        )
        for s in (family.samples if family is not None else [])
    ]
    member_ids = [m.iid for m in members]
    try:
        family_qc = qc_utils.load_family_qc(data_dir, family_id)
    except Exception:
        logger.warning("Failed to read the sample QC of %s", family_id, exc_info=True)
        family_qc = None
    try:
        files = find_ancestry_files(data_dir, family_id)
    except OSError:
        files = None
    try:
        diagnostics = load_family_diagnostics(data_dir, family_id, member_ids)
    except Exception:
        logger.warning("Failed to read the diagnostics of %s", family_id, exc_info=True)
        diagnostics = []

    roh_section, froh = _roh(data_dir, members, family_id)
    sections = [
        _pedigree(members, family_qc),
        _ancestry(files, members, family_id),
        _sample_qc(data_dir, family_id, family_qc, cohort),
        roh_section,
        _pgs(files, members, family_id),
        _aneuploidy(data_dir, members, family_id, family_qc),
        _upd(data_dir, family_id, members, froh, family_qc, diagnostics),
        _diagnostics(data_dir, family_id, diagnostics),
    ]
    return FamilyReport(
        family_id=family_id,
        cohort_name=cohort.name if cohort is not None else "",
        generated_at=datetime.now().isoformat(sep=" ", timespec="minutes"),
        generated_by=user,
        app_version=_app_version(),
        pipeline=_pipeline_text(cohort),
        summary=latest_family_summary(data_dir, family_id),
        sections=sections,
        methods=_methods(sections, files, family_qc, data_dir, family_id, cohort),
    )
