"""What a diagnostic records about its variant when it is saved.

A diagnostic used to store only its key, genes, impact and verdict. The family
report needs more for each finding (HGVS on the MANE Select transcript,
zygosity, segregation, ClinVar, gnomAD, the gene's mode of inheritance, an
SV's copy number and gene overlap), and the table row the curator saves from
already holds it. So the row is copied into :data:`SNAPSHOT_COLUMNS`, after
the diagnostic columns, and the report reads it back without re-reading or
re-filtering the Wombat and SV tables.

The builders are pure: they read a table row and return ``{column: text}``,
``""`` where the row has nothing. ``Source`` is always set, so a record with a
snapshot (even an empty one) is told apart from a record saved before this
existed. STR and UPD diagnostics have none: the report words them from their
own files (:mod:`diagnostic_labels`).
"""

from typing import Any, Dict, Iterable, List, Mapping, Optional
from urllib.parse import unquote

from .inheritance import normalize_gt, observe
from .pedigree_labels import normalize_pheno

#: Columns appended to the diagnostic header, in file order.
SNAPSHOT_COLUMNS: List[str] = [
    "Source", "Transcript", "MANE", "HGVSc", "HGVSp", "Zygosity", "Segregation",
    "Carriers", "ClinVar", "gnomAD_AF", "gnomAD_FAF95", "MOI", "GeneOverlap",
    "CopyNumber",
]  # fmt: skip

#: Pedigree phenotype code -> the word written in ``Carriers``.
_PHENO_WORDS = {"1": "unaffected", "2": "affected"}

#: gnomAD columns of a Wombat table, joint (exomes + genomes) first.
_GNOMAD_AF = ("AF_joint", "AF_genomes", "AF_exomes")
_GNOMAD_FAF = ("fafmax_faf95_max_joint", "fafmax_faf95_max_genomes")

_EMPTY = ("", ".", "-", "None", "nan")


def empty_snapshot() -> Dict[str, str]:
    return dict.fromkeys(SNAPSHOT_COLUMNS, "")


def has_snapshot(record: Mapping[str, Any]) -> bool:
    """Whether a diagnostic record was saved with a snapshot."""
    return bool((record.get("Source") or "").strip())


def _text(value: Any) -> str:
    text = "" if value is None else str(value).strip()
    return "" if text in _EMPTY else text


def _single(value: Any) -> str:
    """One value of an aggregated cell (``a,b`` → ``""``: which one is unknown)."""
    parts = list(dict.fromkeys(p for p in _text(value).split(",") if _text(p)))
    return parts[0] if len(parts) == 1 else ""


def _strip_id(hgvs: str) -> str:
    """``ENST…:c.1A>G`` → ``c.1A>G``; VEP writes ``=`` as ``%3D``."""
    return unquote(hgvs.split(":", 1)[1] if ":" in hgvs else hgvs)


def _pheno_word(code: Optional[str]) -> str:
    return _PHENO_WORDS.get(normalize_pheno(code), "phenotype unknown")


def _carriers(
    key: str,
    sample: str,
    rows: Iterable[Mapping[str, Any]],
    key_of,
    phenotypes: Mapping[str, Optional[str]],
    gt_of,
    extra: Mapping[str, str] = {},
) -> str:
    """Other members carrying ``key``, from the rows given: ``S1 (affected, 0/1)``.

    ``extra`` adds members known from the row itself (the parents' GTs).
    An empty ``phenotypes`` (no pedigree at hand) restricts nothing.
    """
    found: Dict[str, str] = {}
    for other in rows:
        who = str(other.get("sample") or "")
        if not who or who == sample or who in found:
            continue
        if phenotypes and who not in phenotypes:
            continue
        gt = normalize_gt(gt_of(other))
        if key_of(other) == key and gt and any(a not in "0." for a in gt.split("/")):
            found[who] = gt
    for who, gt in extra.items():
        if who and who != sample and who not in found:
            found[who] = gt
    return "; ".join(
        f"{who} ({_pheno_word(phenotypes.get(who))}, {gt})"
        for who, gt in sorted(found.items())
    )


def snv_snapshot(
    row: Mapping[str, Any],
    *,
    source: str,
    rows: Iterable[Mapping[str, Any]] = (),
    phenotypes: Mapping[str, Optional[str]] = {},
    sex: Optional[str] = None,
) -> Dict[str, str]:
    """A Wombat row's snapshot (aggregated by :func:`wombat.aggregate_by_variant`).

    ``rows`` are the other rows of the table, searched for other members
    carrying the variant; ``phenotypes`` maps each family member to its
    pedigree phenotype, and only members are listed. ``sex`` is the sample's,
    when the row has no ``sex`` column.
    """
    snap = empty_snapshot()
    snap["Source"] = source
    # The transcript picked before aggregation; else a row with one transcript.
    hgvsc = _text(row.get("_hgvsc"))
    if hgvsc:
        snap["Transcript"] = _text(row.get("_hgvs_transcript"))
        snap["MANE"] = _text(row.get("_hgvs_mane"))
        hgvsp = _text(row.get("_hgvsp"))
    else:
        hgvsc = _single(row.get("VEP_HGVSc"))
        hgvsp = _single(row.get("VEP_HGVSp")) if hgvsc else ""
        if hgvsc:
            snap["Transcript"] = hgvsc.split(":", 1)[0] if ":" in hgvsc else ""
            snap["MANE"] = _single(row.get("VEP_MANE_SELECT"))
    snap["HGVSc"] = _strip_id(hgvsc) if hgvsc else ""
    snap["HGVSp"] = _strip_id(hgvsp) if hgvsp else ""

    sample = str(row.get("sample") or "")
    father, mother = _text(row.get("father")), _text(row.get("mother"))
    obs = observe(
        chrom=row.get("#CHROM"),
        pos=row.get("POS"),
        sample_gt=row.get("sample_gt"),
        father_gt=row.get("father_gt"),
        mother_gt=row.get("mother_gt"),
        sex=row.get("sex") or sex,
        father_pheno=phenotypes.get(father),
        mother_pheno=phenotypes.get(mother),
    )
    if obs.zygosity != "unknown":
        snap["Zygosity"] = obs.zygosity
    if obs.origin != "unknown":
        snap["Segregation"] = obs.origin

    parents = {}
    for who, gt in ((father, row.get("father_gt")), (mother, row.get("mother_gt"))):
        gt = normalize_gt(gt)
        if who and gt and any(a not in "0." for a in gt.split("/")):
            parents[who] = gt
    key = _snv_key(row)
    snap["Carriers"] = _carriers(
        key, sample, rows, _snv_key, phenotypes,
        lambda r: r.get("sample_gt"), parents,
    )  # fmt: skip

    terms = [
        t.strip().replace("_", " ")
        for part in _text(row.get("VEP_CLIN_SIG")).split(",")
        for t in part.split("&")
        if _text(t)
    ]
    snap["ClinVar"] = ", ".join(dict.fromkeys(terms))
    gnomad = [c for c in (*_GNOMAD_AF, *_GNOMAD_FAF) if c in row]
    af = next((_text(row.get(c)) for c in _GNOMAD_AF if _text(row.get(c))), "")
    faf = next((_text(row.get(c)) for c in _GNOMAD_FAF if _text(row.get(c))), "")
    # A table with gnomAD columns but no value: the site is not in gnomAD.
    snap["gnomAD_AF"] = af or ("absent" if gnomad and not faf else "")
    snap["gnomAD_FAF95"] = faf or ("absent" if gnomad and not af else "")
    snap["MOI"] = moi_text(row)
    return snap


def moi_text(row: Mapping[str, Any]) -> str:
    """The MOI verdicts in words: ``SHANK3: consistent with autosomal dominant``.

    The table's ``MOI`` cell is a badge label (``AD ✓``); the first line of
    each badge's tooltip (:func:`inheritance.annotate_moi`) says it in words.
    """
    lines = [
        str(b.get("tooltip") or "").split("\n")[0].strip()
        for b in row.get("MoiBadges") or []
    ]
    return "; ".join(t for t in lines if t) or _text(row.get("MOI"))


def _snv_key(row: Mapping[str, Any]) -> str:
    return (
        f"{row.get('#CHROM', '')}:{row.get('POS', '')}"
        f":{row.get('REF', '')}:{row.get('ALT', '')}"
    )


def _float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def wisecondorx_snapshot(
    row: Mapping[str, Any], *, source: str = "WisecondorX"
) -> Dict[str, str]:
    """A WisecondorX call's snapshot: genes with their overlap, copy number.

    WisecondorX gives a log2 ratio against its reference; on an autosome the
    copy number is about ``2 × 2^ratio``. On X and Y the expected ploidy
    depends on the karyotype, so only the ratio is written.
    """
    snap = empty_snapshot()
    snap["Source"] = source
    genes = [g.strip() for g in _text(row.get("gene")).split(",") if _text(g)]
    snap["GeneOverlap"] = ", ".join(dict.fromkeys(genes))
    ratio, z = _float(row.get("ratio")), _float(row.get("zscore"))
    if ratio is not None:
        locus = str(row.get("_original_locus") or row.get("chr:start-end") or "")
        chrom = locus.split(":")[0].lower().removeprefix("chr")
        evidence = f"log2 ratio {ratio:+.2f}" + (
            f", z {z:+.1f}" if z is not None else ""
        )
        snap["CopyNumber"] = (
            f"≈{2 * 2**ratio:.1f} ({evidence})"
            if chrom not in ("x", "y", "")
            else evidence
        )
    return snap


def vedette_snapshot(
    row: Mapping[str, Any],
    *,
    source: str,
    rows: Iterable[Mapping[str, Any]] = (),
    phenotypes: Mapping[str, Optional[str]] = {},
) -> Dict[str, str]:
    """A Vedette carrier row's snapshot (:func:`vedette.carrier_rows`)."""
    snap = empty_snapshot()
    snap["Source"] = source
    gt = normalize_gt(row.get("GT"))
    alleles = gt.split("/") if gt else []
    if len(alleles) == 1 and alleles[0] not in "0.":
        snap["Zygosity"] = "hemi"
    elif alleles and all(a not in "0." for a in alleles):
        snap["Zygosity"] = "hom"
    elif any(a not in "0." for a in alleles):
        snap["Zygosity"] = "het"
    snap["Segregation"] = _text(row.get("inheritance"))
    snap["Carriers"] = _carriers(
        str(row.get("site_id") or ""), str(row.get("sample") or ""), rows,
        lambda r: str(r.get("site_id") or ""), phenotypes, lambda r: r.get("GT"),
    )  # fmt: skip
    snap["GeneOverlap"] = _text(row.get("genes"))
    cn = row.get("CN")
    snap["CopyNumber"] = "" if cn is None else str(cn)
    return snap


def snapshot_values(snapshot: Optional[Mapping[str, Any]]) -> List[str]:
    """The snapshot as fields in :data:`SNAPSHOT_COLUMNS` order (all empty for none)."""
    snapshot = snapshot or {}
    return [_text(snapshot.get(c)) for c in SNAPSHOT_COLUMNS]
