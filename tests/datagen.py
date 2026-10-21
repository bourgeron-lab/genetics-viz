"""Build a small synthetic data directory and reference bundle for the tests.

CI has no access to the real data mount, so everything the app reads is
generated here in the formats the pipelines write. The contents are chosen to
exercise the traps found so far:

- family IDs that are numeric for most rows and not for one (``AU-X-1``), the
  case that made polars type ``family_id`` as an integer;
- a pedigree with a header and one without, and float-stringified codes;
  per-family pedigrees are headerless, as the pipeline writes them;
- a cohort directory without a pedigree and one whose pedigree is unreadable;
- a sharded ``families/`` tree;
- two filter-tag variants of one bundle, ``dp6gq15`` and ``dp10gq20``, each
  with admixture proportions (``Q.tsv``);
- a trio whose child has elevated ROH and whose parents are related, so the
  role-dependent ROH wording has something to hide;
- ExpansionHunter tables with a child's expansion above both parents, an AR
  carrier, a male's single chrX allele, a ``LowDepth`` call and a locus with no
  thresholds; family 1002 has no STR results at all;
- a GeneTrek snapshot cache with a renamed gene (``KIAA1109`` -> ``BLTP1``), an
  alias shared by two genes, versioned lists and a gene whose sources disagree
  on its mode of inheritance; and Wombat rows for each MOI verdict;
- CNV coverage bigWigs for the 1001 trio, the carrier halved over its
  WisecondorX loss, a control panel (inside DEMO itself) with a masked stretch
  and PAR-aware sex groups, reached through a compute-cluster
  ``cnv_panel_dir``; and a reference-free CRAM for S3 with one of each read
  signature, split reads and a soft-clip cluster at the loss's breakpoints;
- sample QC (somalier) tables for 1001 and AU-X-1 and the DEMO cohort: S1 is
  not S3's father, S3 is contaminated, S5 reads XXY, S6 is male in the data
  but female in the pedigree, and S2 is related to S5 across families; 1002
  has none. ``validations/sampleqc.tsv`` acknowledges two of them: S6's sex
  (explained) and the S2-S5 relatedness, whose note names a relationship; an
  earlier acknowledgement of S3's contamination was withdrawn.
"""

import gzip
import hashlib
import json
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from genetics_viz.utils.sharding import compute_shard_prefix

BUNDLE = "1.1.0"
TAGS = ("dp6gq15", "dp10gq20")
REGIONS = ("EUR", "CSA")

#: (FID, IID, father, mother, sex, phenotype)
DEMO_PEDIGREE: List[Sequence[str]] = [
    ("1001", "S1", "0", "0", "1", "1"),
    ("1001", "S2", "0", "0", "2", "1"),
    ("1001", "S3", "S1", "S2", "1", "2"),
    ("1002", "S4", "0", "0", "2.0", "-9.0"),
    ("AU-X-1", "S5", "-9", "", "1", "2"),
    ("AU-X-1", "S6", "S5", "0", "2", "1"),
]

FROH_COLUMNS = [
    "IID", "region", "af_source", "call_rate", "qc_pass", "nroh", "sroh_mb",
    "nroh_1.5mb", "largest_roh_mb", "largest_roh_cm", "froh_0.5mb", "froh_1.5mb",
    "froh_5mb", "froh_10mb", "f_short", "f_intermediate", "f_long", "ref_group",
    "ref_n", "pct_froh_1.5mb", "pct_f_short", "pct_f_intermediate", "pct_f_long",
    "ref_median_froh_1.5mb", "homozygosity_pattern", "froh_category",
    "upd_candidate_chrom", "father", "mother", "parents_kinship", "expected_f",
]  # fmt: skip

METRICS = [
    "nroh", "sroh_mb", "nroh_1.5mb", "largest_roh_mb", "largest_roh_cm",
    "froh_0.5mb", "froh_1.5mb", "froh_5mb", "froh_10mb", "f_short",
    "f_intermediate", "f_long",
]  # fmt: skip

TRAITS = ("height_test", "bmi_test")


def write_tsv(
    path: Path, header: Sequence[str], rows: Iterable[Sequence], gz: bool = False
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = (
        "\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows)
    )
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def family_dir(data: Path, fid: str) -> Path:
    s1, s2 = compute_shard_prefix(fid)
    return data / "families" / s1 / s2 / fid


def sample_dir(data: Path, sid: str) -> Path:
    s1, s2 = compute_shard_prefix(sid)
    return data / "samples" / s1 / s2 / sid


# ---------------------------------------------------------------------------
# Reference bundle
# ---------------------------------------------------------------------------


def build_reference(root: Path) -> Path:
    """``<root>/v1.1.0/bundle`` with PCA labels and the ROH tier."""
    bundle = root / f"v{BUNDLE}" / "bundle"
    pcs_rows, froh_rows = [], []
    for r, region in enumerate(REGIONS):
        for i in range(30):
            iid = f"REF_{region}_{i}"
            pcs_rows.append(
                [iid, round(r * 0.05 + i * 0.001, 5), round(-r * 0.03 + i * 0.0005, 5),
                 0.001 * i, -0.001 * i, f"POP{region}", region]
            )  # fmt: skip
            short = 0.03 + 0.0005 * i + (0.02 if region == "CSA" else 0)
            inter = 0.002 + 0.0002 * i
            long = 0.0 if i < 25 else 0.004 * (i - 24)
            froh = inter + long
            froh_rows.append(
                [iid, f"POP{region}", region, region, 150, 120.5, 9, 4.2, 6.1,
                 short + froh, froh, long, long, short, inter, long]
            )  # fmt: skip
    write_tsv(
        bundle / "labels" / "reference_pcs.tsv.gz",
        ["IID", "PC1", "PC2", "PC3", "PC4", "population", "region"],
        pcs_rows,
        gz=True,
    )
    write_tsv(
        bundle / "roh" / "reference_froh.tsv.gz",
        ["IID", "population", "region", "af_source", *METRICS],
        froh_rows,
        gz=True,
    )
    write_tsv(
        bundle / "roh" / "hotspots.tsv.gz",
        ["group", "chrom", "start", "end", "max_freq"],
        [("ALL", "1", 45300000, 46200000, 0.36), ("EUR", "2", 1000000, 1300000, 0.3)],
        gz=True,
    )
    return root


# ---------------------------------------------------------------------------
# Per-family results
# ---------------------------------------------------------------------------


def _froh_row(iid: str, father: str, mother: str, elevated: bool) -> List:
    if elevated:
        values = dict(froh=0.105, short=0.035, inter=0.012, long=0.093, cat="second-degree-like",
                      pattern="recent", pct=96.5, kin=0.07)  # fmt: skip
    else:
        values = dict(froh=0.008, short=0.04, inter=0.006, long=0.0, cat="none",
                      pattern="typical", pct=45.0, kin=None)  # fmt: skip
    kin = "NA" if values["kin"] is None else values["kin"]
    return [
        iid, "EUR", "EUR", 0.999, "true", 140, 120.0, 8, 73.8 if elevated else 4.0,
        40.1 if elevated else 2.0, values["short"] + values["froh"], values["froh"],
        values["long"], values["long"], values["short"], values["inter"], values["long"],
        "EUR", 30, values["pct"], 50.0, 40.0, 96.0 if elevated else 44.0, 0.007,
        values["pattern"], values["cat"], ".", father, mother, kin, kin,
    ]  # fmt: skip


def build_family_results(
    data: Path, fid: str, members: Sequence[Sequence[str]]
) -> None:
    """The family pedigree, ancestry (two tags) and ROH (one tag) results."""
    fdir = family_dir(data, fid)
    fdir.mkdir(parents=True, exist_ok=True)
    (fdir / f"{fid}.pedigree.tsv").write_text(
        "".join("\t".join(m) + "\n" for m in members)
    )
    for tag in TAGS:
        prefix = fdir / "ancestry" / f"{fid}.apgs_b{BUNDLE}_{tag}"
        write_tsv(
            Path(f"{prefix}.pcs.tsv"),
            ["IID", *[f"PC{i}" for i in range(1, 21)]],
            [[m[1], *[0.01 * (k + 1) for k in range(20)]] for m in members],
        )
        write_tsv(
            Path(f"{prefix}.ancestry.tsv"),
            ["IID", "region", "population", "region_confidence", "population_confidence",
             "mean_knn_distance", "nearest_reference_distance", "local_reference_density",
             "is_outlier"],
            [[m[1], "EUR", "POPEUR", 0.98, 0.9, 0.01, 0.005, 12.0, "False"] for m in members],
        )  # fmt: skip
        write_tsv(
            Path(f"{prefix}.pgs_zscore.tsv"),
            ["IID", *TRAITS],
            [[m[1], 0.5, "nan"] for m in members],
        )
        write_tsv(
            Path(f"{prefix}.Q.tsv"),
            ["IID", "snp_coverage", *REGIONS],
            [[m[1], 0.99, 0.8, 0.2] for m in members],
        )

    rprefix = fdir / "roh" / f"{fid}.apgs_b{BUNDLE}_{TAGS[0]}"
    children = {m[1] for m in members if m[2] not in ("0", "-9", "")}
    write_tsv(
        Path(f"{rprefix}.froh.tsv"),
        FROH_COLUMNS,
        [_froh_row(m[1], m[2] if m[2] not in ("0", "") else ".", m[3] if m[3] not in ("0", "") else ".",
                   elevated=m[1] in children) for m in members],
    )  # fmt: skip
    write_tsv(
        Path(f"{rprefix}.roh.tsv"),
        ["IID", "chrom", "start", "end", "length_mb", "length_cm", "n_sites", "quality",
         "roh_class", "hotspot_group", "hotspot_fraction", "hotspot_fraction_pooled"],
        [[m[1], "chr2", 1000001, 9000000, 8.0, 12.0, 4000, 90, "long", "EUR", 0.0, 0.0]
         for m in members],
    )  # fmt: skip
    write_tsv(
        Path(f"{rprefix}.roh_genes.tsv"),
        ["IID", "chrom", "roh_start", "roh_end", "roh_length_mb", "roh_length_cm", "roh_class",
         "gene", "gene_id", "gene_start", "gene_end", "gene_fraction_in_roh",
         "gene_hotspot_fraction", "hotspot_group"],
        [[m[1], "chr2", 1000001, 9000000, 8.0, 12.0, "long", "GENEA", "ENSG0001", 2000000,
          2100000, 1.0, 0.0, "EUR"] for m in members],
    )  # fmt: skip
    if len(members) > 1:
        ids = [m[1] for m in members]
        pairs = [(a, b) for i, a in enumerate(ids) for b in ids[i + 1 :]]
        parents = {(m[2], m[3]) for m in members}

        def phi(a: str, b: str) -> float:
            return 0.07 if (a, b) in parents or (b, a) in parents else 0.25

        write_tsv(
            Path(f"{rprefix}.kinship.tsv"),
            ["IID1", "IID2", "n_sites", "ibs0", "kinship"],
            [[a, b, 1040000, 0.0001, phi(a, b)] for a, b in pairs],
        )


# ---------------------------------------------------------------------------
# Cohorts
# ---------------------------------------------------------------------------


def build_cohort_results(
    cohort_dir: Path, name: str, pedigree: Sequence[Sequence[str]]
) -> None:
    """Cohort-level ancestry and ROH tables, with the extra ``family_id`` column."""
    for tag in TAGS:
        prefix = cohort_dir / "ancestry" / f"{name}.apgs_b{BUNDLE}_{tag}"
        write_tsv(
            Path(f"{prefix}.pcs.tsv"),
            ["IID", *[f"PC{i}" for i in range(1, 21)], "family_id"],
            [[m[1], *[0.01 * (k + 1) for k in range(20)], m[0]] for m in pedigree],
        )
        write_tsv(
            Path(f"{prefix}.ancestry.tsv"),
            ["IID", "region", "population", "region_confidence", "population_confidence",
             "mean_knn_distance", "nearest_reference_distance", "local_reference_density",
             "is_outlier", "family_id"],
            [[m[1], "EUR", "POPEUR", 0.98, 0.9, 0.01, 0.005, 12.0, "False", m[0]] for m in pedigree],
        )  # fmt: skip
        write_tsv(
            Path(f"{prefix}.pgs_zscore.tsv"),
            ["IID", *TRAITS, "family_id"],
            [[m[1], 0.1 * i, "nan", m[0]] for i, m in enumerate(pedigree)],
        )
        write_tsv(
            Path(f"{prefix}.Q.tsv"),
            ["IID", "snp_coverage", *REGIONS, "family_id"],
            [[m[1], 0.99, 0.8, 0.2, m[0]] for m in pedigree],
        )
    rprefix = cohort_dir / "roh" / f"{name}.apgs_b{BUNDLE}_{TAGS[0]}"
    children = {m[1] for m in pedigree if m[2] not in ("0", "-9", "")}
    write_tsv(
        Path(f"{rprefix}.froh.tsv"),
        [*FROH_COLUMNS, "family_id"],
        [
            [*_froh_row(m[1], ".", ".", elevated=m[1] in children), m[0]]
            for m in pedigree
        ],
    )
    write_tsv(
        Path(f"{rprefix}.roh_windows.tsv"),
        ["group", "chrom", "start", "end", "n_founders_in_roh", "n_founders", "freq"],
        [("ALL", "chr1", 45300000, 45400000, 2, 4, 0.5), ("ALL", "chr1", 45400000, 45500000, 1, 4, 0.25),
         ("EUR", "chr2", 1000000, 1100000, 1, 4, 0.25)],
    )  # fmt: skip


def build_state(
    cohort_dir: Path, name: str, pedigree_file: Path, params_file: Path
) -> None:
    def sha(p: Path) -> str:
        return hashlib.sha256(p.read_bytes()).hexdigest()

    record = {
        "status": "success",
        "run_id": "run-1",
        "started_at": "2026-10-01T10:00:00+00:00",
        "finished_at": "2026-10-01T12:30:00+00:00",
        "cohort_name": name,
        "pipeline": {"version": "1.2.0", "revision": "main", "commit_id": "abcdef1234567"},
        "pedigree": {"path": "/pasteur/x.tsv", "sha256": sha(pedigree_file), "families": 3, "individuals": 6},
        "params_file": {"path": "/pasteur/x.yml", "sha256": sha(params_file)},
        "steps_requested": ["ancestry", "roh"],
        "completion_measured": "before_run",
        "completion": {"ancestry": {"done": 3, "total": 3, "pct": 100.0}, "roh": None},
    }  # fmt: skip
    (cohort_dir / ".ghfc-ngs.state.json").write_text(
        json.dumps({"schema_version": 1, "cohort_name": name, "last_run": record,
                    "last_successful_run": record, "history": [record]})
    )  # fmt: skip


# ---------------------------------------------------------------------------
# Short tandem repeats (ExpansionHunter)
# ---------------------------------------------------------------------------

EH_LABEL = "jmsa50"

#: eh_family_table's COLUMNS, in order.
EH_COLUMNS = [
    "family_id", "sample", "father", "mother", "sex", "phenotype",
    "locus", "chrom", "pos", "end", "motif", "ref_units", "filter", "genotype",
    "allele1", "allele2", "allele1_ci", "allele2_ci",
    "allele1_support", "allele2_support",
    "allele1_spanning", "allele2_spanning", "allele1_flanking", "allele2_flanking",
    "allele1_inrepeat", "allele2_inrepeat", "locus_coverage",
    "disease", "inheritance", "normal_max", "intermediate_min", "pathogenic_min",
    "allele1_class", "allele2_class", "status", "parents_max_allele",
    "premutation_min", "report", "read_length", "allele1_confidence",
    "allele2_confidence", "severity",
]  # fmt: skip

#: locus -> (chrom, pos, end, motif, ref_units, disease, inheritance,
#: normal_max, intermediate_min, premutation_min, pathogenic_min, report);
#: "." where the catalog has none.
EH_LOCI = {
    "AFF3": ("chr2", 100104799, 100104822, "GCC", 8, ".", ".", ".", ".", ".", ".", "."),
    "HTT": ("chr4", 3074876, 3074933, "CAG", 19, "Huntington disease", "AD", 26, 27,
            ".", 36, "flag"),
    "FXN": ("chr9", 69037286, 69037304, "GAA", 6, "Friedreich ataxia", "AR", 30, ".",
            ".", 70, "flag"),
    "TCF4": ("chr18", 55586155, 55586227, "CTG", 25, "Fuchs dystrophy", "AD", 40, ".",
             ".", 50, "risk"),
    "ARX_2": ("chrX", 25013529, 25013565, "NGC", 12, "DEE1", "XLR", 12, ".", ".", 20,
              "flag"),
    "FMR1": ("chrX", 147912050, 147912110, "CGG", 20, "Fragile X", "XLD", 40, 41, 55,
             201, "flag"),
}  # fmt: skip

#: (sample, locus) -> (allele1, allele2, support1, support2) where it differs
#: from the default: 10 and 12 units, spanning (one allele on a male's chrX).
#: HTT in S3 is pathogenic (high), FMR1 in S2 a premutation (medium), FXN in
#: S2 a carrier and TCF4 in S5 a risk allele (info). ARX_2 in S6 is pathogenic
#: but of low confidence at a GCN motif with no spanning read (medium).
EH_CALLS = {
    ("S3", "HTT"): (18, 40, "SPANNING", "SPANNING"),
    ("S2", "FXN"): (8, 120, "SPANNING", "INREPEAT"),
    ("S2", "FMR1"): (30, 60, "SPANNING", "INREPEAT"),
    ("S5", "TCF4"): (14, 60, "SPANNING", "INREPEAT"),
    ("S6", "AFF3"): (9, 31, "SPANNING", "FLANKING"),
    ("S6", "ARX_2"): (12, 22, "SPANNING", "FLANKING"),
}

#: Read length per family, from the members' ExpansionHunter JSON: AU-X-1 is
#: an older 100 bp genome.
EH_READ_LENGTH = {"1001": 150, "AU-X-1": 100}

_SEX = {"1": "male", "2": "female"}
_PHENO = {"1": "unaffected", "2": "affected"}


def _eh_class(count, normal, inter, premut, patho) -> str:
    if normal == "." and inter == "." and premut == "." and patho == ".":
        return "."
    if patho != "." and count >= patho:
        return "pathogenic"
    if premut != "." and count >= premut:
        return "premutation"
    if inter != "." and count >= inter:
        return "intermediate"
    if normal == "." or count <= normal:
        return "normal"
    return "above_normal"


def _eh_confidence(count, cls, spanning, inrepeat, motif, ref, patho, length) -> str:
    """eh_family_table's confidence of one allele past normal."""
    if count == "." or cls in (".", "normal"):
        return "."
    span = length - 20
    if "N" in motif or count * len(motif) <= span:
        return "ok" if spanning >= 2 else "low"
    if ref * len(motif) > span and patho != "." and patho * len(motif) < 2 * length:
        return "low"
    return "ok" if inrepeat >= 3 else "low"


def _eh_severity(status, confidences, report) -> str:
    """eh_family_table's severity of a row."""
    rank = {"pathogenic": 0, "premutation": 1, "intermediate": 2, "carrier": 2}.get(
        status
    )
    if rank is None:
        return "."
    if report == "risk":
        rank = 2
    if "low" in confidences:
        rank += 1
    return ("high", "medium", "info")[min(rank, 2)]


def _eh_rows(fid: str, members: Sequence[Sequence[str]]) -> List[List]:
    """One row per locus and member, as eh_family_table writes them."""
    ids = {m[1] for m in members}
    rows: List[List] = []
    longest: Dict = {}
    length = EH_READ_LENGTH.get(fid, 150)
    for locus, spec in EH_LOCI.items():
        chrom, pos, end, motif, ref, disease, inh, n, i, pm, p, report = spec
        for m in members:
            sid, sex = m[1], _SEX.get(m[4].split(".")[0], ".")
            a1, a2, s1, s2 = EH_CALLS.get(
                (sid, locus), (10, 12, "SPANNING", "SPANNING")
            )
            single = chrom == "chrX" and sex == "male"
            if single:
                a1, a2, s2 = a2, ".", "."
            classes = [_eh_class(a, n, i, pm, p) if a != "." else "." for a in (a1, a2)]
            ranked = [c for c in classes if c != "."]
            if not ranked:
                status = "."
            elif inh == "AR" and len(ranked) == 2 and ranked.count("pathogenic") == 1:
                status = "carrier"
            else:
                order = ["normal", "above_normal", "intermediate", "premutation",
                         "pathogenic"]  # fmt: skip
                status = max(ranked, key=order.index)
            father = m[2] if m[2] in ids else "."
            mother = m[3] if m[3] in ids else "."
            parents = [
                longest[(x, locus)] for x in (father, mother) if (x, locus) in longest
            ]
            low = sid == "S6" and locus == "ARX_2"
            # (spanning, in-repeat) reads of each allele
            reads = [(12, 0)]
            if a2 != ".":
                reads.append(
                    (9 if s2 == "SPANNING" else 0, 4 if s2 == "INREPEAT" else 0)
                )
            reads.append((0, 0))
            confidences = [
                _eh_confidence(a, c, r[0], r[1], motif, ref, p, length)
                if inh != "."
                else "."
                for a, c, r in zip((a1, a2), classes, reads)
            ]
            severity = _eh_severity(status, confidences, report) if inh != "." else "."
            rows.append([
                fid, sid, father, mother, sex, _PHENO.get(m[5].split(".")[0], "."),
                locus, chrom, pos, end, motif, ref, "LowDepth" if low else "PASS",
                "1" if single else "0/1",
                a1, a2, f"{a1}-{a1}", "." if a2 == "." else f"{a2}-{a2 + 3}", s1, s2,
                reads[0][0], 0 if a2 == "." else reads[1][0], 30,
                0 if a2 == "." else 28, 0, 0 if a2 == "." else reads[1][1],
                8.5 if low else 31.2, disease, inh, n, i, p, classes[0], classes[1],
                status, max(parents) if parents else ".", pm, report, length,
                confidences[0], confidences[1], severity,
            ])  # fmt: skip
            longest[(sid, locus)] = max(a for a in (a1, a2) if a != ".")
    return rows


def build_str_results(
    data: Path, fid: str, members: Sequence[Sequence[str]]
) -> List[List]:
    """A family's ExpansionHunter table and one REViewer image. Returns the rows."""
    rows = _eh_rows(fid, members)
    write_tsv(
        family_dir(data, fid)
        / "str"
        / "expansionhunter"
        / f"{fid}.{EH_LABEL}.genotypes.tsv",
        EH_COLUMNS,
        rows,
    )
    for m in members:
        sdir = sample_dir(data, m[1]) / "str" / "expansionhunter"
        write_tsv(sdir / f"{m[1]}.{EH_LABEL}.vcf", ["#CHROM"], [])
        # EXPANSIONHUNTER_REVIEWER: a metrics table for every member, an image
        # for flagged genotypes only. S2's flagged FXN is left undrawn.
        rdir = sdir / "reviewer"
        write_tsv(
            rdir / f"{m[1]}.{EH_LABEL}.metrics.tsv",
            ["VariantId", "Genotype", "AlleleDepth"],
            [],
        )
        if m[1] == "S3":
            (rdir / f"S3.{EH_LABEL}.HTT.svg").write_text(
                '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"/>'
            )
    return rows


def build_cohort_str_results(cohort_dir: Path, name: str, rows: List[List]) -> None:
    """eh_cohort_merge's three tables, from the families' rows."""
    col = {c: k for k, c in enumerate(EH_COLUMNS)}

    def longest(r):
        return max(a for a in (r[col["allele1"]], r[col["allele2"]]) if a != ".")

    founders = [r for r in rows if r[col["father"]] == "." and r[col["mother"]] == "."]
    by_locus: Dict[str, List[int]] = {}
    for r in founders:
        by_locus.setdefault(r[col["locus"]], []).append(longest(r))

    def pct(r):
        values = by_locus[r[col["locus"]]]
        return round(100 * sum(v <= longest(r) for v in values) / len(values), 1)

    edir = cohort_dir / "str" / "expansionhunter"
    write_tsv(
        edir / f"{name}.{EH_LABEL}.genotypes.tsv",
        [*EH_COLUMNS, "cohort_percentile"],
        [[*r, pct(r)] for r in rows],
    )
    loci_rows = []
    statuses = ("pathogenic", "premutation", "carrier", "intermediate", "above_normal")
    severities = ("high", "medium", "info")
    for locus, spec in EH_LOCI.items():
        chrom, pos, end, motif, _ref, disease, inh, n, i, pm, p, report = spec
        values = sorted(by_locus[locus])
        mine = [r for r in rows if r[col["locus"]] == locus]
        count = {s: sum(r[col["status"]] == s for r in mine) for s in statuses}
        count.update(
            {s: sum(r[col["severity"]] == s for r in mine) for s in severities}
        )
        loci_rows.append([
            locus, chrom, pos, end, motif, disease, inh, n, i, pm, p, report, len(mine),
            len(values), values[len(values) // 2], values[-1], values[-1], values[-1],
            values[-1], *(count[s] for s in statuses + severities),
            sum(r[col["filter"]] != "PASS" for r in mine),
        ])  # fmt: skip
    write_tsv(
        edir / f"{name}.{EH_LABEL}.loci.tsv",
        ["locus", "chrom", "pos", "end", "motif", "disease", "inheritance", "normal_max",
         "intermediate_min", "premutation_min", "pathogenic_min", "report",
         "n_samples", "n_founders", "founders_median", "founders_p95", "founders_p99",
         "founders_p99.9", "founders_max", *(f"n_{s}" for s in statuses + severities),
         "n_not_pass"],
        loci_rows,
    )  # fmt: skip
    write_tsv(
        edir / f"{name}.{EH_LABEL}.flags.tsv",
        [*EH_COLUMNS, "cohort_percentile", "founders_median", "founders_p99.9",
         "founders_max"],
        [[*r, pct(r), 10, 12, 12] for r in rows
         if r[col["severity"]] in ("high", "medium")],
    )  # fmt: skip


# WisecondorX's cohort aberrations: one robust loss, one robust gain.
WCX_COLUMNS = [
    "chr", "start", "end", "ratio", "zscore", "type", "source", "barcode",
    "genic_symbol", "genic_ensg", "exonic_symbol", "exonic_ensg", "family_id",
]  # fmt: skip
WCX_ROWS = [
    ["chr2", 165114001, 165316000, -1.0329, -9.6, "primary", "loss", "S3", "SCN2A",
     "ENSG00000136531", "SCN2A", "ENSG00000136531", "1001"],
    ["chr2", 200000001, 200500000, 0.4939, 6.0, "primary", "gain", "S6", "", "", "",
     "", "AU-X-1"],
]  # fmt: skip


def build_data_dir(root: Path) -> Path:
    """The whole synthetic data directory. Returns its path."""
    data = root / "data"
    cohorts = data / "cohorts"

    demo = cohorts / "DEMO"
    pedigree_file = write_tsv(
        demo / "DEMO.pedigree.tsv",
        ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"],
        DEMO_PEDIGREE,
    )
    params_file = demo / "DEMO.params.yml"
    params_file.write_text(
        "# test params\ncohort_name: DEMO\nsteps:\n  - ancestry\n  - roh\n"
        f"cnv_panel_dir: {CNV_PANEL_CLUSTER_PATH}\n"
    )
    build_state(demo, "DEMO", pedigree_file, params_file)
    build_cohort_results(demo, "DEMO", DEMO_PEDIGREE)
    write_tsv(
        demo / "svs" / "wisecondorx" / "DEMO_aberrations.bed", WCX_COLUMNS, WCX_ROWS
    )

    # Same data without a header line.
    nohead = cohorts / "NOHEAD"
    nohead.mkdir(parents=True)
    (nohead / "NOHEAD.pedigree.tsv").write_text(
        "".join("\t".join(r) + "\n" for r in DEMO_PEDIGREE[:3])
    )

    (cohorts / "EMPTY").mkdir()  # no pedigree at all
    broken = cohorts / "BROKEN"
    broken.mkdir()
    (broken / "BROKEN.pedigree.tsv").write_text("only_one_column\nvalue\n")

    by_family: Dict[str, List[Sequence[str]]] = {}
    for row in DEMO_PEDIGREE:
        by_family.setdefault(row[0], []).append(row)
    str_rows: List[List] = []
    for fid, members in by_family.items():
        build_family_results(data, fid, members)
        if fid != "1002":
            str_rows += build_str_results(data, fid, members)
    build_cohort_str_results(demo, "DEMO", str_rows)
    # The family's own WisecondorX calls, for its SVs tab.
    write_tsv(
        family_dir(data, "1001")
        / "svs"
        / "wisecondorx"
        / "1001_aberrations.annotated.bed",
        WCX_COLUMNS,
        [r for r in WCX_ROWS if r[-1] == "1001"],
    )
    build_cohort_qc_results(
        demo,
        "DEMO",
        [build_qc_results(data, fid, by_family[fid]) for fid in QC_FAMILIES],
    )
    build_qc_acks(data)
    build_wombat_results(data)
    build_cnv_results(data, demo)
    return data


# ---------------------------------------------------------------------------
# GeneTrek snapshot and Wombat results
# ---------------------------------------------------------------------------

GENETREK_RELEASE = "test-release"

#: GeneTrek's field list, trimmed to what the app reads, in the API's shape.
GENETREK_FIELDS: List[Dict] = [
    {"name": "hgnc_symbol", "category": "Gene identity", "value_type": "string", "id_type": "symbol"},
    {"name": "hgnc_name", "category": "Gene identity", "value_type": "string"},
    {"name": "hgnc_prev_symbol", "category": "Gene identity", "value_type": "string"},
    {"name": "hgnc_alias_symbol", "category": "Gene identity", "value_type": "string"},
    {"name": "ensg", "category": "Gene identity", "value_type": "string", "id_type": "ensg"},
    {"name": "loeuf", "category": "Constraint", "value_type": "float"},
    {"name": "eagle_score", "category": "Gene lists", "value_type": "categorical"},
    {"name": "hc_ndd_v7", "fullname": "High-confidence NDD v7", "category": "Gene lists",
     "subcategory": "GeneTrek", "visibility": "primary", "value_type": "boolean",
     "aliases": ["hc_ndd"], "n_genes": 3, "description": "High-confidence NDD genes"},
    {"name": "hc_ndd_v6", "fullname": "High-confidence NDD v6", "category": "Gene lists",
     "subcategory": "GeneTrek", "visibility": "secondary", "value_type": "boolean",
     "n_genes": 2},
    {"name": "sfari_1", "fullname": "SFARI 1", "category": "Gene lists",
     "subcategory": "Autism", "visibility": "secondary", "value_type": "boolean",
     "aliases": [], "n_genes": 1},
    {"name": "constraint", "fullname": "Constrained genes", "category": "Constraint",
     "visibility": "primary", "value_type": "boolean", "n_genes": 2,
     "derived_from": ["loeuf"]},
    {"name": "syngo", "fullname": "SynGO", "category": "Function",
     "visibility": "primary", "value_type": "boolean", "n_genes": 1},
    {"name": "inheritance", "fullname": "Mode of inheritance", "category": "Inheritance",
     "value_type": "categorical"},
    {"name": "inheritance_panelapp_t", "fullname": "Mode of inheritance from PanelApp",
     "category": "Inheritance", "value_type": "categorical"},
]  # fmt: skip

#: symbol, name, prev, alias, ensg, loeuf, eagle, hc7, hc6, sfari, constraint,
#: syngo, inheritance, inheritance_panelapp
GENETREK_GENES: List[Sequence] = [
    ("SHANK3", "SH3 and ankyrin 3", None, "PROSAP2", "ENSG00000251322", 0.36, "Definitive",
     True, True, True, True, True, "Autosomal Dominant", "Autosomal or pseudoautosomal dominant"),
    ("BLTP1", "bridge-like lipid transfer 1", "KIAA1109", None, "ENSG00000138688", 0.62, None,
     True, True, False, False, False, "Autosomal Recessive", None),
    ("XRG", "an X-linked recessive gene", None, None, "ENSG00000000001", None, None,
     True, False, False, False, False, "X-linked Recessive", None),
    ("ARONLY", "a recessive gene", None, "SHARED", "ENSG00000000002", 0.9, None,
     False, False, False, True, False, "Autosomal Recessive", None),
    ("BOTH", "dominant in one source", None, "SHARED", "ENSG00000000003", 0.9, None,
     False, False, False, False, False, "Autosomal Recessive", "Autosomal or pseudoautosomal dominant"),
    ("PLAIN", "in no set", None, None, "ENSG00000000004", None, None,
     False, False, False, False, False, None, None),
]  # fmt: skip


def genetrek_frame():
    """GeneTrek's table for the genes above, typed as the API's parquet is."""
    import polars as pl

    names = [
        "hgnc_symbol", "hgnc_name", "hgnc_prev_symbol", "hgnc_alias_symbol", "ensg",
        "loeuf", "eagle_score", "hc_ndd_v7", "hc_ndd_v6", "sfari_1", "constraint",
        "syngo", "inheritance", "inheritance_panelapp_t",
    ]  # fmt: skip
    df = pl.DataFrame([list(g) for g in GENETREK_GENES], schema=names, orient="row")
    return df.with_columns(
        pl.col("loeuf").cast(pl.Float32),
        pl.col("eagle_score", "inheritance", "inheritance_panelapp_t").cast(
            pl.Categorical
        ),
    )


def build_genetrek_cache(cache_dir: Path) -> Path:
    """The on-disk cache the app writes after a download."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    genetrek_frame().write_parquet(cache_dir / f"{GENETREK_RELEASE}.parquet")
    (cache_dir / "fields.json").write_text(
        json.dumps(
            {
                "release": GENETREK_RELEASE,
                "fetched_at": "2026-10-01T12:00:00",
                "fields": GENETREK_FIELDS,
            }
        )
    )
    return cache_dir


WOMBAT_COLUMNS = [
    "#CHROM", "POS", "REF", "ALT", "sample", "sex", "sample_gt", "father",
    "mother", "father_gt", "mother_gt", "VEP_Consequence", "VEP_SYMBOL",
    "VEP_Gene", "fafmax_faf95_max_genomes",
]  # fmt: skip

#: Family 1001's child S3 (male, affected; parents S1, S2 unaffected):
#: an AD de novo, an AR compound het in trans (one row under the gene's old
#: symbol), an X-linked recessive hemizygous call, an AR single het, and a
#: gene GeneTrek gives no mode.
WOMBAT_ROWS = [
    ("chr22", 50674000, "G", "A", "S3", "1", "0/1", "S1", "S2", "0/0", "0/0",
     "missense_variant", "SHANK3", "ENSG00000251322", 0.0),
    ("chr7", 32000000, "C", "T", "S3", "1", "0/1", "S1", "S2", "0/1", "0/0",
     "stop_gained", "KIAA1109", "ENSG00000138688", 0.0001),
    ("chr7", 32100000, "A", "G", "S3", "1", "0|1", "S1", "S2", "0/0", "0/1",
     "missense_variant", "BLTP1", "ENSG00000138688", 0.0002),
    ("chrX", 70000000, "T", "C", "S3", "1", "1/1", "S1", "S2", "0/0", "0/1",
     "missense_variant", "XRG", "ENSG00000000001", 0.0),
    ("chr1", 1000, "A", "C", "S3", "1", "0/1", "S1", "S2", "0/1", "0/0",
     "synonymous_variant", "ARONLY", "ENSG00000000002", 0.001),
    ("chr2", 2000, "G", "T", "S3", "1", "0/1", "S1", "S2", "0/0", "0/1",
     "missense_variant", "PLAIN", "ENSG00000000004", 0.001),
]  # fmt: skip


def build_wombat_results(data: Path) -> None:
    """Family 1001's Wombat tables (two configs) and the DEMO cohort's two."""
    wombat = family_dir(data, "1001") / "wombat"
    write_tsv(
        wombat / "1001.rare.vep1.annotated.default.tsv", WOMBAT_COLUMNS, WOMBAT_ROWS
    )
    # A second config, so the family page has a second Wombat sub-tab
    write_tsv(
        wombat / "1001.rare.vep1.annotated.strict.tsv",
        WOMBAT_COLUMNS,
        WOMBAT_ROWS[:2],
    )
    # A de novo config: the sample QC of 1001 puts a warning above it. Sorted
    # after "default", so the first sub-tab is unchanged.
    write_tsv(
        wombat / "1001.rare.vep1.annotated.rare_denovo.tsv",
        WOMBAT_COLUMNS,
        WOMBAT_ROWS[:1],
    )
    write_tsv(
        data / "cohorts" / "DEMO" / "wombat" / "DEMO.rare.vep1.default.results.tsv",
        WOMBAT_COLUMNS,
        WOMBAT_ROWS,
    )
    write_tsv(
        data / "params" / "genesets" / "local_list.tsv",
        ["gene"],
        [("ARONLY",), ("PLAIN",)],
    )
    build_headerless_cohort_wombat(data)


def build_headerless_cohort_wombat(data: Path) -> Path:
    """DEMO's "strict" cohort table as several real cohorts' de novo tables are.

    The older ".annotated." name, and no header: a blank line, then the rows.
    Family 1001's "strict" table lends it its header. Part of the default data;
    returns the path for tests that rewrite it or its family file.
    """
    path = (
        data
        / "cohorts"
        / "DEMO"
        / "wombat"
        / "DEMO.rare.vep1.annotated.strict.results.tsv"
    )
    path.write_text(
        "\n" + "".join("\t".join(map(str, r)) + "\n" for r in WOMBAT_ROWS[:2])
    )
    return path


# ---------------------------------------------------------------------------
# CNV coverage, control panel and read evidence
# ---------------------------------------------------------------------------

CNV_TAG = "cov_b250q1"
CNV_BIN = 250
#: As a cohort's params file writes it: a path on the compute cluster.
CNV_PANEL_CLUSTER_PATH = "/pasteur/helix/projects/x/WGS/cohorts/DEMO/cnv_panel/v1"
CNV_CHROMS = {"chr2": 242193529, "chrX": 156040895, "chrY": 57227415}
#: Only these stretches are written: the whole grid would be millions of bins.
CNV_SPANS = {
    "chr2": (164_800_000, 165_700_000),
    "chrX": (2_700_000, 2_900_000),
    "chrY": (2_800_000, 2_900_000),
}
#: S3's WisecondorX loss (WCX_ROWS), and a stretch where the panel is unusable.
CNV_LOSS = (165_114_001, 165_316_000)
CNV_MASKED = (165_400_000, 165_410_000)
PAR1_END = 2_781_479
CNV_SEX = {"S1": "XY", "S2": "XX", "S3": "XY"}
#: About the fixture CRAM's depth outside the loss, so its read depth reads 2.
CNV_MEDIAN_DEPTH = 9.0
#: The reference-free CRAM decodes without it; the file only has to exist.
CRAM_REFERENCE = "cram_reference/GRCh38.fasta"


def _bigwig(path: Path, value_at) -> None:
    import pybigtools

    def entries():
        for chrom, (lo, hi) in CNV_SPANS.items():
            for start in range(lo, hi, CNV_BIN):
                value = value_at(chrom, start)
                if value is not None:
                    yield chrom, start, start + CNV_BIN, float(value)

    path.parent.mkdir(parents=True, exist_ok=True)
    pybigtools.open(str(path), "w").write(CNV_CHROMS, entries())


def _sample_value(sid: str):
    sex = CNV_SEX[sid]

    def value(chrom: str, start: int):
        if chrom == "chr2":
            in_loss = sid == "S3" and CNV_LOSS[0] <= start < CNV_LOSS[1]
            return 1.0 if in_loss else 2.0
        if chrom == "chrX":
            return 2.0 if sex == "XX" or start < PAR1_END else 1.0
        return 1.0 if sex == "XY" else 0.0

    return value


def _panel_value(group: str, stat: str):
    def value(chrom: str, start: int):
        if group == "autosomes" and chrom != "chr2":
            return None
        if group != "autosomes" and chrom == "chr2":
            return None
        if group == "XX" and chrom != "chrX":
            return None
        if stat == "mad":
            return 0.1
        if chrom == "chr2":
            return 0.4 if CNV_MASKED[0] <= start < CNV_MASKED[1] else 2.0
        if chrom == "chrX":
            return 2.0 if group == "XX" or start < PAR1_END else 1.0
        return 1.0

    return value


def cram_reference(root: Path) -> Path:
    return root / CRAM_REFERENCE


def build_cnv_results(data: Path, demo: Path) -> None:
    for sid in CNV_SEX:
        cov = sample_dir(data, sid) / "coverage"
        _bigwig(cov / f"{sid}.{CNV_TAG}.bw", _sample_value(sid))
        (cov / f"{sid}.{CNV_TAG}.json").write_text(
            json.dumps(
                {
                    "sample": sid,
                    "tag": CNV_TAG,
                    "bin": CNV_BIN,
                    "mapq": 1,
                    "autosomal_median_depth": CNV_MEDIAN_DEPTH,
                }
            )
        )

    panel = demo / "cnv_panel" / "v1"
    for group in ("autosomes", "XX", "XY"):
        for stat in ("median", "mad"):
            _bigwig(
                panel / f"panel.{CNV_TAG}.{group}.{stat}.bw", _panel_value(group, stat)
            )
    par = {
        "build": "GRCh38",
        "coordinates": "1-based inclusive",
        "chrX": {"PAR1": [10001, PAR1_END], "PAR2": [155701383, 156030895]},
        "chrY": {"PAR1": [10001, PAR1_END], "PAR2": [56887903, 57217415]},
    }
    (panel / "panel.json").write_text(
        json.dumps(
            {
                "panel_version": "1",
                "tag": CNV_TAG,
                "bin": CNV_BIN,
                "mapq": 1,
                "reference": "GRCh38.fasta",
                "source_cohort": "DEMO",
                "n_controls": {"autosomes": 6, "XX": 3, "XY": 3},
                "mad_scale": 1.0,
                "max_fit_mad": 0.3,
                "par": par,
                "expected_ploidy": {
                    "autosomes": 2,
                    "XX": {"chrX": 2},
                    "XY": {"chrX": {"non_par": 1, "par": 2}, "chrY": {"non_par": 1}},
                },
            }
        )
    )
    write_tsv(
        demo / "coverage" / f"DEMO.{CNV_TAG}.panel_v1.fit.tsv",
        ["family_id", "sample", "pedigree_sex", "coverage_sex", "fit_mad",
         "chrX_ratio", "chrY_ratio"],
        [["1001", "S1", "male", "XY", 0.12, 0.5, 0.5],
         ["1001", "S2", "female", "XX", 0.13, 1.0, 0.0],
         ["1001", "S3", "male", "XY", 0.41, 0.5, 0.5]],
    )  # fmt: skip

    ref = cram_reference(data.parent)
    ref.parent.mkdir(parents=True, exist_ok=True)
    ref.write_text(">chr2\nACGT\n")
    # Lengths as in the CRAM header, so htslib does not warn about them.
    ref.with_suffix(".fasta.fai").write_text(
        "".join(f"{c}\t{n}\t6\t4\t5\n" for c, n in CNV_CHROMS.items())
    )
    build_cram(sample_dir(data, "S3") / "sequences" / "S3.GRCh38_GIABv3.cram")


#: Signatures written at S3's loss: (name, kind, left start, right start).
CRAM_SIGNATURES = (
    # Fragments that fit the insert size across the junction (about 330 bp).
    [(f"del{i}", "deletion", CNV_LOSS[0] - 230 + 20 * i, CNV_LOSS[1] + 10 * i)
     for i in range(6)]
    + [(f"inv{i}", "inversion_pp", CNV_LOSS[0] - 600 + 30 * i, CNV_LOSS[1] - 300)
       for i in range(2)]
    + [("dup0", "duplication", CNV_LOSS[0] - 300, CNV_LOSS[0] + 100)]
)  # fmt: skip
CRAM_SPLITS = 2
CRAM_CLIPPED = 3
#: Concordant pairs every 22 bp over +-1.7 kb of each breakpoint, every other
#: one dropped inside the (heterozygous) loss: depth about 9x, halved inside.
CRAM_CONCORDANT_STEP = 22


def concordant_starts(side: int) -> List[int]:
    centre = CNV_LOSS[side]
    starts = []
    for i, pos in enumerate(range(centre - 1700, centre + 1700, CRAM_CONCORDANT_STEP)):
        inside = pos >= centre if side == 0 else pos + 350 <= centre
        if inside and i % 2:
            continue
        starts.append(pos)
    return starts


CRAM_CONCORDANT_LEFT = len(concordant_starts(0))
CRAM_CONCORDANT_PER_SIDE = 60


def build_cram(path: Path) -> None:
    """A reference-free CRAM: pairs of every class at both breakpoints."""
    import pysam

    read_len = 100
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": n} for c, n in CNV_CHROMS.items()]
        + [{"SN": "chr5", "LN": 181538259}],
    }
    chrom_id = {c: i for i, c in enumerate([*CNV_CHROMS, "chr5"])}
    reads: List = []

    def add(name, flag, chrom, pos, mate_chrom, mate_pos, tlen, cigar="100M",
            tags=(), mapq=60):  # fmt: skip
        reads.append((chrom_id[chrom], pos, name, flag, mate_chrom, mate_pos, tlen,
                      cigar, tags, mapq))  # fmt: skip

    def pair(name, pos, mate_pos, left_flag, right_flag, mapq=60):
        tlen = mate_pos + read_len - pos
        mq = (("MQ", mapq),)
        add(name, left_flag, "chr2", pos, "chr2", mate_pos, tlen, tags=mq, mapq=mapq)
        add(name, right_flag, "chr2", mate_pos, "chr2", pos, -tlen, tags=mq,
            mapq=mapq)  # fmt: skip

    for side, centre in enumerate(CNV_LOSS):
        for i, pos in enumerate(concordant_starts(side)):
            pair(f"c{side}_{i}", pos, pos + 250 + (i % 5) * 10, 99, 147)
    for name, kind, left, right in CRAM_SIGNATURES:
        flags = {
            "deletion": (97, 145),  # +-, not proper
            "inversion_pp": (65, 129),  # ++
            "duplication": (81, 161),  # -+
        }[kind]
        pair(name, left, right, *flags)
    # A deletion pair in a repeat: MAPQ 0, shown only when asked for.
    pair("lowdel0", CNV_LOSS[0] - 500, CNV_LOSS[1] + 300, 97, 145, mapq=0)
    add("trans0", 97, "chr2", CNV_LOSS[0] - 200, "chr5", 1_000_000, 0)
    # Split and clipped reads all end at the left breakpoint: one clip cluster.
    for i in range(CRAM_CLIPPED):
        add(f"clip{i}", 0, "chr2", CNV_LOSS[0] - 70, "chr2", CNV_LOSS[0] - 70, 0,
            cigar="70M30S")  # fmt: skip
    for i in range(CRAM_SPLITS):
        add(f"split{i}", 0, "chr2", CNV_LOSS[0] - 60, "chr2",
            CNV_LOSS[0] - 300, 0, cigar="60M40S",
            tags=(("SA", f"chr2,{CNV_LOSS[1] + 1},+,60S40M,60,0;"),))  # fmt: skip

    path.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(
        str(path), "wc", header=header, format_options=[b"no_ref=1"]
    ) as out:
        for (
            rid,
            pos,
            name,
            flag,
            mate_chrom,
            mate_pos,
            tlen,
            cigar,
            tags,
            mapq,
        ) in sorted(reads, key=lambda r: (r[0], r[1])):
            seg = pysam.AlignedSegment(out.header)
            seg.query_name = name
            seg.query_sequence = "A" * read_len
            seg.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            seg.flag = flag
            seg.reference_id = rid
            seg.reference_start = pos
            seg.mapping_quality = mapq
            seg.cigarstring = cigar
            seg.next_reference_id = chrom_id[mate_chrom]
            seg.next_reference_start = mate_pos
            seg.template_length = tlen
            for tag, value in tags:
                seg.set_tag(tag, value)
            out.write(seg)
    pysam.index(str(path))


# ---------------------------------------------------------------------------
# Sample QC (somalier)
# ---------------------------------------------------------------------------

QC_LABEL = "somalier_hg38"
QC_FAMILIES = ("1001", "AU-X-1")

QC_SAMPLE_COLUMNS = [
    "family_id", "sample", "father", "mother", "pedigree_sex", "inferred_sex",
    "karyotype", "sex_evidence", "x_copies", "y_copies", "somalier_x_copies",
    "somalier_y_copies", "x_het_rate", "x_sites", "depth_mean", "n_sites_called",
    "het_hom_alt_ratio", "contamination_charr", "charr_sites", "somalier_charr",
    "p_middling_ab", "trio_informative_sites", "trio_mendelian_error_rate", "flags",
]  # fmt: skip
QC_PAIR_COLUMNS = [
    "family_id_a", "sample_a", "family_id_b", "sample_b", "n_sites", "ibs0",
    "ibs2", "ibs0_ibs2_ratio", "king_kinship", "somalier_relatedness",
    "hom_concordance", "x_ibs0", "x_ibs2", "expected_kinship",
    "expected_relationship", "observed_relationship", "po_informative_sites",
    "po_opposite_hom_rate", "status",
]  # fmt: skip
QC_FLAG_COLUMNS = [
    "level", "family_id", "sample_a", "sample_b", "flag", "severity", "evidence",
]  # fmt: skip

#: sample -> (inferred sex, karyotype, X, Y, X het, CHARR, trio sites, trio rate)
QC_SAMPLES = {
    "S1": ("male", "XY", 0.99, 0.91, 0.0, 0.012, None, None),
    "S2": ("female", "XX", 1.99, 0.0, 0.62, 0.018, None, None),
    "S3": ("male", "XY", 0.99, 0.90, 0.0, 0.031, 4500000, 0.0008),
    "S5": ("male", "XXY", 1.98, 0.90, 0.35, 0.015, None, None),
    "S6": ("male", "XY", 1.00, 0.92, 0.0, 0.019, None, None),
}

#: (sample_a, sample_b, kinship, ibs0/ibs2, expected kinship, expected,
#: observed, PO opposite-hom rate, status)
QC_PAIRS = {
    "1001": [
        ("S1", "S2", 0.0600, 0.1400, 0.0, "unrelated", "third-degree", None,
         "mismatch"),
        ("S1", "S3", 0.0051, 0.1700, 0.25, "parent-offspring", "unrelated", 0.0300,
         "mismatch"),
        ("S2", "S3", 0.2512, 0.0, 0.25, "parent-offspring", "parent-offspring",
         0.0004, "ok"),
    ],
    "AU-X-1": [
        ("S5", "S6", 0.2490, 0.0, 0.25, "parent-offspring", "parent-offspring", None,
         "ok"),
    ],
}  # fmt: skip

QC_FLAGS = {
    "1001": [
        ("S1", "S3", "NOT_PARENT", "high",
         "S1 is not the father of S3: expected parent-offspring (kinship 0.2500),"
         " observed unrelated (kinship 0.0051, ibs0/ibs2 0.1700, 16200 sites);"
         " opposite homozygotes at 0.0300 of 4000000 sites in the family BCF"),
        ("S1", "S2", "CONSANGUINITY", "info",
         "spouses: expected unrelated (kinship 0.0000), observed third-degree"
         " (kinship 0.0600, ibs0/ibs2 0.1400, 16200 sites)"),
        ("S3", ".", "POSSIBLE_CONTAMINATION", "medium",
         "CHARR 0.031 at 8200 sites: about 3.1% of the reads from another person"),
    ],
    "AU-X-1": [
        ("S5", ".", "KARYOTYPE_XXY", "high", "copy numbers X 1.98 Y 0.90 (coverage)"),
        ("S6", ".", "SEX_MISMATCH", "high",
         "pedigree female, inferred male (XY, copy numbers X 1.00 Y 0.92 (coverage))"),
    ],
}  # fmt: skip

#: The cross-family pairs the cohort sees: S2 (1001) and S5 (AU-X-1), second
#: degree with the IBS0 of relatives (high), and S1 (1001) and S5 again,
#: third degree (information).
QC_CROSS_PAIRS = [
    ("1001", "S2", "AU-X-1", "S5", 0.1100, 0.0500, "second-degree", "high", ""),
    ("1001", "S1", "AU-X-1", "S5", 0.0500, 0.1400, "third-degree", "info",
     "; the third degree is within the reach of shared ancestry or admixture"),
]  # fmt: skip

_PED_SEX = {"1": "male", "2": "female"}


def _qc_parent(parent: str) -> str:
    """The pipeline writes ``0`` for every unknown-parent spelling."""
    return "0" if parent in ("", "0", "-9", ".") else parent


def _f4(value) -> str:
    """The pipeline's float format: four decimals, ``.`` for missing."""
    return "." if value is None else f"{float(value):.4f}"


def _qc_pair_row(fa, a, fb, b, phi, ratio, exp_phi, expected, observed, po, status):
    ibs2 = 9800
    return [
        fa, a, fb, b, 16200, int(round(ratio * ibs2)), ibs2, _f4(ratio), _f4(phi),
        _f4(phi * 2), ".", 0, 200, _f4(exp_phi), expected, observed,
        4000000 if po is not None else ".", _f4(po), status,
    ]  # fmt: skip


def build_qc_results(
    data: Path, fid: str, members: Sequence[Sequence[str]]
) -> Dict[str, List[List]]:
    """``<family>/qc/<FID>.somalier_hg38.*``. Returns the rows for the cohort."""
    qc = family_dir(data, fid) / "qc"
    flags = [
        ["family", fid, a, b, code, sev, ev] for a, b, code, sev, ev in QC_FLAGS[fid]
    ]
    by_sample: Dict[str, List[str]] = {}
    for row in flags:
        for sample in (row[2], row[3]):
            if sample != ".":
                by_sample.setdefault(sample, []).append(row[4])
    samples = []
    for _, sid, father, mother, sex, _pheno in members:
        inferred, karyotype, x, y, het, charr, trio_n, trio_rate = QC_SAMPLES[sid]
        samples.append(
            [
                fid, sid, _qc_parent(father), _qc_parent(mother),
                _PED_SEX.get(sex.replace(".0", ""), "unknown"), inferred, karyotype,
                "coverage", _f4(x), _f4(y), _f4(x), _f4(y), _f4(het), 330, "35.0000",
                16800, "1.9000", _f4(charr), 8200, _f4(4 * charr), "0.0230",
                trio_n if trio_n is not None else ".", _f4(trio_rate),
                ",".join(sorted(set(by_sample.get(sid, [])))) or ".",
            ]
        )  # fmt: skip
    pairs = [_qc_pair_row(fid, a, fid, b, *rest) for a, b, *rest in QC_PAIRS[fid]]
    write_tsv(qc / f"{fid}.{QC_LABEL}.samples.tsv", QC_SAMPLE_COLUMNS, samples)
    write_tsv(qc / f"{fid}.{QC_LABEL}.pairs.tsv", QC_PAIR_COLUMNS, pairs)
    write_tsv(qc / f"{fid}.{QC_LABEL}.flags.tsv", QC_FLAG_COLUMNS, flags)
    (qc / f"{fid}.{QC_LABEL}.inputs.tsv").write_text(
        f"#inputs\tv1\n#step\tsampleqc\n#label\t{QC_LABEL}\n#family\t{fid}\n"
    )
    return {"samples": samples, "pairs": pairs, "flags": flags}


def build_cohort_qc_results(
    cohort_dir: Path, name: str, families: List[Dict[str, List[List]]]
) -> None:
    """The cohort tables: the families', plus the cross-family pairs and their flags."""
    qc = cohort_dir / "qc"
    cross_flags = [
        ["cohort", f"{fa},{fb}", a, b, "CRYPTIC_RELATED", severity,
         f"families {fa} and {fb}: observed {observed} (kinship {phi:.4f},"
         f" ibs0/ibs2 {ratio:.4f}, 16200 sites){note}"]
        for fa, a, fb, b, phi, ratio, observed, severity, note in QC_CROSS_PAIRS
    ]  # fmt: skip
    related = {s for p in QC_CROSS_PAIRS for s in (p[1], p[3])}
    samples = [list(r) for fam in families for r in fam["samples"]]
    for row in samples:
        if row[1] in related:
            codes = set(row[-1].split(",")) - {"."} | {"CRYPTIC_RELATED"}
            row[-1] = ",".join(sorted(codes))
    pairs = [r for fam in families for r in fam["pairs"]] + [
        _qc_pair_row(fa, a, fb, b, phi, ratio, 0.0, "unrelated", observed, None,
                     "mismatch")
        for fa, a, fb, b, phi, ratio, observed, _sev, _note in QC_CROSS_PAIRS
    ]  # fmt: skip
    flags = [r for fam in families for r in fam["flags"]] + cross_flags
    write_tsv(qc / f"{name}.{QC_LABEL}.samples.tsv", QC_SAMPLE_COLUMNS, samples)
    write_tsv(qc / f"{name}.{QC_LABEL}.pairs.tsv", QC_PAIR_COLUMNS, pairs)
    write_tsv(qc / f"{name}.{QC_LABEL}.flags.tsv", QC_FLAG_COLUMNS, flags)
    (qc / f"{name}.{QC_LABEL}.html").write_text("<html></html>")


QC_ACK_COLUMNS = [
    "Flag", "FID", "SampleA", "SampleB", "Status", "Note", "User", "Timestamp",
]  # fmt: skip

#: The rows of ``validations/sampleqc.tsv``: the last one of a key wins.
QC_ACKS = [
    ["POSSIBLE_CONTAMINATION", "1001", "S3", "", "confirmed", "", "curator",
     "2026-09-01T10:00:00"],
    ["SEX_MISMATCH", "AU-X-1", "S6", "", "explained",
     "Pedigree entry error, fixed at the next release", "curator",
     "2026-09-02T10:00:00"],
    ["CRYPTIC_RELATED", "1001,AU-X-1", "S2", "S5", "explained",
     "Known first cousins: their mothers are sisters", "curator",
     "2026-09-03T10:00:00"],
    ["POSSIBLE_CONTAMINATION", "1001", "S3", "", "withdrawn", "", "curator",
     "2026-09-04T10:00:00"],
    ["CONSANGUINITY", "1001", "S1", "S2", "explained",
     "Declared: the parents are second cousins", "curator", "2026-09-05T10:00:00"],
]  # fmt: skip


def build_qc_acks(data: Path) -> None:
    """Curator acknowledgements of sample QC flags, as the app appends them."""
    write_tsv(data / "validations" / "sampleqc.tsv", QC_ACK_COLUMNS, QC_ACKS)


# ---------------------------------------------------------------------------
# Uniparental disomy (opt-in: it would change every flag count of 1001)
# ---------------------------------------------------------------------------

#: A UPD_CANDIDATE row as ``sampleqc_classify`` words it: S2 (the mother) and
#: S3 are opposite homozygotes on chr2, so S3 has two paternal copies.
UPD_FLAG = (
    "S2", "S3", "UPD_CANDIDATE", "high",
    "96 of the 100 opposite homozygotes of S2 and S3 lie on chr2, between 0.1 and"
    " 39.7 Mb; S3 is heterozygous at 0.004 of its 812 sites there (0.310"
    " genome-wide), at 1.00 of its mean depth: uniparental isodisomy: both copies"
    " from S3's other parent",
)  # fmt: skip


def add_upd_case(data: Path, fid: str = "1001", chrom: str = "2") -> None:
    """Give 1001's S3 a UPD of ``chrom``: the ROH ``upd_candidate_chrom`` and
    the Sample QC flag, in the family files and the cohort's."""
    a, b, code, severity, evidence = UPD_FLAG
    evidence = evidence.replace("chr2", f"chr{chrom}")
    for path in [
        *(family_dir(data, fid) / "qc").glob(f"{fid}.{QC_LABEL}.flags.tsv"),
        *(data / "cohorts").glob(f"*/qc/*.{QC_LABEL}.flags.tsv"),
    ]:
        if path.name.startswith(f"{fid}.") or f"\t{fid}\t" in path.read_text():
            with open(path, "a") as f:
                f.write(
                    "\t".join(["family", fid, a, b, code, severity, evidence]) + "\n"
                )
    for path in [
        *(family_dir(data, fid) / "roh").glob("*.froh.tsv"),
        *(data / "cohorts").glob("*/roh/*.froh.tsv"),
    ]:
        lines = path.read_text().splitlines()
        header = lines[0].split("\t")
        col, iid = header.index("upd_candidate_chrom"), header.index("IID")
        fam = header.index("family_id") if "family_id" in header else None
        for i, line in enumerate(lines[1:], 1):
            cells = line.split("\t")
            if cells[iid] == b and (fam is None or cells[fam] == fid):
                cells[col] = chrom
                lines[i] = "\t".join(cells)
        path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Family report
# ---------------------------------------------------------------------------

#: The polygenic scores the family report shows (``config/report.yaml``).
REPORT_TRAITS = (
    "autism_grove2026", "adhd_demontis2023", "epilepsy_ILAE2023/all", "int_savage2018",
)  # fmt: skip

#: z-scores per sample, in ``REPORT_TRAITS`` order; S3's ADHD score is missing.
REPORT_Z = {
    "S1": (0.2, -0.4, 0.0, 1.1),
    "S2": (-0.6, 0.3, 0.5, 0.4),
    "S3": (2.4, "nan", 0.5, -1.2),
}

DIAGNOSTIC_COLUMNS = [
    "FID", "Variant", "Gene", "Impact", "Sample", "User", "Timestamp", "Comment",
    "Ignore", "Diagnostic",
]  # fmt: skip

#: (file, row): a pathogenic SNV and an uncertain SV in S3, a benign SNV and an
#: ignored pathogenic one, which the report leaves out.
REPORT_DIAGNOSTICS = [
    ("snvs.tsv", ["1001", "chr22:50721500:G:A", "SHANK3", "stop_gained", "S3",
                  "curator", "2026-10-01T10:00:00", "de novo, fits the phenotype",
                  "0", "pathogenic"]),
    ("svs.tsv", ["1001", "chr2:165114001-165316000:del", "SCN2A:exonic", "LOSS", "S3",
                 "curator", "2026-10-02T10:00:00", "", "0", "uncertain"]),
    ("snvs.tsv", ["1001", "chr1:1000:C:T", "PLAIN", "missense_variant", "S1",
                  "curator", "2026-10-03T10:00:00", "", "0", "benign"]),
    ("snvs.tsv", ["1001", "chr1:2000:C:T", "PLAIN", "missense_variant", "S2",
                  "curator", "2026-10-04T10:00:00", "", "1", "pathogenic"]),
]  # fmt: skip

VALIDATION_COLUMNS = [
    "FID", "Variant", "Sample", "User", "Inheritance", "Validation", "Comment",
    "Ignore", "Timestamp",
]  # fmt: skip


def _add_columns(path: Path, names: Sequence[str], values: Dict[str, Sequence]) -> None:
    """Append ``names`` to a TSV, the row's values looked up by its first cell."""
    lines = path.read_text().splitlines()
    out = ["\t".join([lines[0], *names])]
    for line in lines[1:]:
        iid = line.split("\t")[0]
        extra = values.get(iid, ("nan",) * len(names))
        out.append("\t".join([line, *(str(v) for v in extra)]))
    path.write_text("\n".join(out) + "\n")


def add_report_case(data: Path, fid: str = "1001") -> None:
    """What the family report reads beyond the default data, for ``fid``.

    - the report's four polygenic scores, in the family and cohort files;
    - a whole-chromosome 21 gain for S3 in the family's WisecondorX calls;
    - diagnostics (see :data:`REPORT_DIAGNOSTICS`) and the SNV's validation.

    Opt-in, as :func:`add_upd_case`: it changes the PGS tables other tests
    compare whole.
    """
    for path in [
        *(family_dir(data, fid) / "ancestry").glob("*.pgs_zscore.tsv"),
        *(data / "cohorts").glob("*/ancestry/*.pgs_zscore.tsv"),
    ]:
        _add_columns(path, REPORT_TRAITS, REPORT_Z)
    wcx = (
        family_dir(data, fid)
        / "svs"
        / "wisecondorx"
        / f"{fid}_aberrations.annotated.bed"
    )
    with open(wcx, "a") as f:
        f.write("\t".join(["chr21", "10500001", "46700000", "0.58", "12.0", "primary",
                           "gain", "S3", "", "", "", "", fid]) + "\n")  # fmt: skip
    for name, row in REPORT_DIAGNOSTICS:
        path = data / "diagnostics" / name
        if not path.exists():
            write_tsv(path, DIAGNOSTIC_COLUMNS, [])
        with open(path, "a") as f:
            f.write("\t".join(row) + "\n")
    write_tsv(
        data / "validations" / "snvs.tsv",
        VALIDATION_COLUMNS,
        [[fid, "chr22:50721500:G:A", "S3", "curator", "de novo", "present", "", "0",
          "2026-10-01T09:00:00"]],
    )  # fmt: skip
