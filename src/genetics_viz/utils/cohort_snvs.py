"""The cohort SNVs/Indels tab without UI: de novo candidates across families.

The ``denovo`` step's cohort ``dnm`` table holds every de novo candidate of
the cohort, one row per variant and child, with its flags, strict verdict,
origin and recurrence. This module turns it into the tab's two tables:

- the candidates of the children shown, filtered by strictness, impact and
  validation status, with the curators' validations and diagnostics;
- the genes hit in two or more families, the cohort-level de novo signal.

Inherited variants are not here: they are the Search page's.
"""

import re
from typing import Any, Collection, Dict, Iterable, List, Optional, Set, Tuple

from . import denovo
from .validation_status import STATUS_COLORS, UNVALIDATED
from .vep import get_consequence_impact, get_highest_consequence_term

Key = Tuple[str, str]

#: VEP impacts, most severe first; the filter's options.
IMPACTS = ("HIGH", "MODERATE", "LOW", "MODIFIER")
#: Shown by default: what a curator reviews first.
DEFAULT_IMPACTS = ("HIGH", "MODERATE")
#: Impact badge colours, from the most severe.
IMPACT_COLORS = {
    "HIGH": "#c62828",
    "MODERATE": "#ef6c00",
    "LOW": "#2e7d32",
    "MODIFIER": "#757575",
}
#: Separators of the dnm table's ``genes`` cell (VEP joins with ``,`` or ``&``).
_GENE_SPLIT = re.compile(r"[,;&|]")


def split_genes(value: Any) -> List[str]:
    """The distinct gene symbols of a ``genes`` cell, in order."""
    if not value:
        return []
    genes = (g.strip() for g in _GENE_SPLIT.split(str(value)))
    return list(dict.fromkeys(g for g in genes if g and g != "."))


def impact_of(consequences: Any) -> str:
    """The most severe VEP impact of a row's consequence terms, else ``""``."""
    terms = (t.strip() for t in re.split(r"[,&]", str(consequences or "")))
    impacts = {get_consequence_impact(t) for t in terms if t and t != "."}
    return next((i for i in IMPACTS if i in impacts), "")


def candidate_rows(
    rows: Iterable[Dict[str, Any]],
    keys: Optional[Collection[Key]],
    curation: Optional[denovo.Curation] = None,
) -> List[Dict[str, Any]]:
    """The candidates of the children in ``keys`` (``None``: all), as table rows.

    Each row keeps the dnm fields the table and the export need, plus
    ``Variant`` (``chr:pos:ref:alt``), ``impact``, ``consequence`` (the most
    severe term), ``_genes``, the badges, and the curators' status (the
    validations' consensus, else ``TODO``) and diagnostic when ``curation`` is
    given. Sorted by family, child and position.
    """
    out: List[Dict[str, Any]] = []
    for r in rows:
        fid, child = str(r.get("family_id") or ""), str(r.get("child") or "")
        if keys is not None and (fid, child) not in keys:
            continue
        variant, _ = denovo.candidate_key(r)
        consequences = r.get("consequences") or ""
        flags = denovo.split_flags(r.get("flags"))
        row: Dict[str, Any] = {
            "_key": f"{fid}:{child}:{variant}",
            "FID": fid,
            "child": child,
            "Variant": variant,
            "var_type": r.get("var_type") or "",
            "genes": r.get("genes") or "",
            "_genes": split_genes(r.get("genes")),
            "consequence": get_highest_consequence_term(consequences)
            if consequences
            else "",
            "impact": impact_of(consequences),
            "strict": str(r.get("strict") or ""),
            "flags": ",".join(flags),
            "flags_badges": denovo.flag_badges(flags),
            "origin": r.get("origin") or "",
            "child_vaf": r.get("child_vaf"),
            "recurrence_families": r.get("recurrence_families"),
            "_pos": _int(r.get("pos")),
            "_chrom": str(r.get("chrom") or ""),
        }
        row["impact_badges"] = (
            [{"label": row["impact"], "color": IMPACT_COLORS[row["impact"]]}]
            if row["impact"]
            else []
        )
        row["strict_badges"] = denovo.strict_badges(row["strict"])
        if curation is not None:
            status, inheritance = curation.status((fid, variant, child))
            status = status or UNVALIDATED
            row["Validation"] = status
            row["Validation_badges"] = [
                {
                    "label": status if not inheritance else f"{status} · {inheritance}",
                    "color": STATUS_COLORS.get(status, "#9e9e9e"),
                }
            ]
            row["Diagnostic"] = curation.diagnostics.get((fid, variant, child), "")
        out.append(row)
    out.sort(key=lambda r: (r["FID"], r["child"], r["_chrom"], r["_pos"]))
    return out


def _int(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def filter_candidates(
    rows: List[Dict[str, Any]],
    *,
    strict_only: bool = True,
    impacts: Collection[str] = DEFAULT_IMPACTS,
    statuses: Collection[str] = (),
) -> List[Dict[str, Any]]:
    """The candidates the tab's filters leave.

    ``strict_only`` keeps the strict ones (``yes``); an empty ``impacts`` or
    ``statuses`` filters nothing, as elsewhere in the app.
    """
    wanted_impacts, wanted_statuses = set(impacts), set(statuses)
    return [
        r
        for r in rows
        if (not strict_only or r["strict"] == "yes")
        and (not wanted_impacts or r["impact"] in wanted_impacts)
        and (not wanted_statuses or r.get("Validation", UNVALIDATED) in wanted_statuses)
    ]


def recurrent_genes(
    rows: Iterable[Dict[str, Any]], minimum_families: int = 2
) -> List[Dict[str, Any]]:
    """Genes hit by candidates of ``minimum_families`` families or more.

    One row per gene: its families and children, how many candidates, their
    consequences and the worst impact. Most families first, then the gene.
    """
    by_gene: Dict[str, Dict[str, Any]] = {}
    for r in rows:
        for gene in r["_genes"]:
            entry = by_gene.setdefault(
                gene,
                {
                    "families": set(),
                    "children": set(),
                    "n": 0,
                    "terms": {},
                    "impacts": set(),
                },
            )
            entry["families"].add(r["FID"])
            entry["children"].add((r["FID"], r["child"]))
            entry["n"] += 1
            entry["terms"][r["consequence"]] = None
            if r["impact"]:
                entry["impacts"].add(r["impact"])
    out = []
    for gene, e in by_gene.items():
        if len(e["families"]) < minimum_families:
            continue
        worst = next((i for i in IMPACTS if i in e["impacts"]), "")
        out.append(
            {
                "gene": gene,
                "n_families": len(e["families"]),
                "n_children": len(e["children"]),
                "n_candidates": e["n"],
                "families": ", ".join(sorted(e["families"])),
                "consequences": ", ".join(t for t in e["terms"] if t),
                "impact": worst,
                "impact_badges": (
                    [{"label": worst, "color": IMPACT_COLORS[worst]}] if worst else []
                ),
            }
        )
    out.sort(key=lambda g: (-g["n_families"], -g["n_candidates"], g["gene"]))
    return out


def families_of(rows: Iterable[Dict[str, Any]]) -> Set[str]:
    return {r["FID"] for r in rows}
