"""Build a small synthetic data directory and reference bundle for the tests.

CI has no access to the real data mount, so everything the app reads is
generated here in the formats the pipelines write. The contents are chosen to
exercise the traps found so far:

- family IDs that are numeric for most rows and not for one (``AU-X-1``), the
  case that made polars type ``family_id`` as an integer;
- a pedigree with a header and one without, and float-stringified codes;
- a cohort directory without a pedigree and one whose pedigree is unreadable;
- a sharded ``families/`` tree;
- two filter-tag variants of one bundle, ``dp6gq15`` and ``dp10gq20``;
- a trio whose child has elevated ROH and whose parents are related, so the
  role-dependent ROH wording has something to hide.
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
    """Ancestry (two tags) and ROH (one tag) results for one family."""
    fdir = family_dir(data, fid)
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
    )
    build_state(demo, "DEMO", pedigree_file, params_file)
    build_cohort_results(demo, "DEMO", DEMO_PEDIGREE)

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
    for fid, members in by_family.items():
        build_family_results(data, fid, members)
    return data
