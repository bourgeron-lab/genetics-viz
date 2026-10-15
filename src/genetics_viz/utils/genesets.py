"""Gene sets for filtering variant tables.

Three sources, offered together in the geneset filter menus (Search page,
family Wombat tab):

- **My tiers**: the genes the current user's GeneTrek tiers colour, e.g.
  "Tiers 1-2" (see ``utils/gene_scoring``);
- **GeneTrek**: every boolean GeneTrek gene list, from the cached snapshot;
- **This data directory**: ``<data_dir>/params/genesets/*.tsv``, one gene
  symbol per line after a header line.

Membership is matched on upper-cased identifiers. GeneTrek sets are expanded
with the unambiguous previous and alias symbols of their genes, so a VEP
``KIAA1109`` matches a set holding ``BLTP1``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, List, Optional, Set

import polars as pl

from genetics_viz.utils.gene_scoring import GeneScorer
from genetics_viz.utils.user_preferences import N_TIERS, SHOW

#: Gene separators inside one cell: VEP joins with ``&``, aggregation with ``,``.
_GENE_TOKEN = r"[^,&\s]+"

GROUP_TIERS = "My tiers"
GROUP_LOCAL = "This data directory"


@dataclass(frozen=True)
class GenesetSource:
    """One selectable gene set."""

    key: str
    label: str
    group: str
    #: Upper-cased identifiers matching the set (aliases included).
    symbols: FrozenSet[str]
    n_genes: int
    description: str = ""
    #: Older versions and sub-lists, shown on demand.
    secondary: bool = False


def load_genesets(data_dir: Path) -> Dict[str, Set[str]]:
    """Load gene sets from params/genesets/*.tsv.

    Args:
        data_dir: Root data directory containing params/genesets/.

    Returns:
        Dictionary mapping geneset name (file stem) to a set of
        upper-cased gene symbols.
    """
    genesets_dir = data_dir / "params" / "genesets"
    result: Dict[str, Set[str]] = {}
    if not genesets_dir.exists():
        return result
    for geneset_file in genesets_dir.glob("*.tsv"):
        genes: Set[str] = set()
        with open(geneset_file, "r") as f:
            next(f, None)  # skip header
            for line in f:
                gene = line.strip()
                if gene:
                    genes.add(gene.upper())
        if genes:
            result[geneset_file.stem] = genes
    return result


def geneset_sources(
    data_dir: Optional[Path], scorer: GeneScorer
) -> List[GenesetSource]:
    """Every gene set a filter menu offers, in display order."""
    snap = scorer.snapshot
    sources: List[GenesetSource] = []

    previous: FrozenSet[str] = frozenset()
    for max_tier in range(1, N_TIERS + 1):
        symbols = scorer.tier_symbols(max_tier)
        if len(symbols) == len(previous):
            continue  # an empty tier would repeat the one above
        previous = symbols
        label = "Tier 1" if max_tier == 1 else f"Tiers 1–{max_tier}"
        sources.append(
            GenesetSource(
                key=f"tier:{max_tier}",
                label=label,
                group=GROUP_TIERS,
                symbols=snap.matching_identifiers(symbols),
                n_genes=len(symbols),
                description="Genes your profile colours at this tier or better",
            )
        )

    # The user's own sets first within each category, then GeneTrek's order.
    def used(name: str) -> bool:
        a = scorer.assignments.get(name)
        return isinstance(a, int) or a == SHOW

    for gs in sorted(snap.catalog, key=lambda g: not used(g.name)):
        members = snap.set_members.get(gs.name, frozenset())
        sources.append(
            GenesetSource(
                key=f"gt:{gs.key}",
                label=gs.label,
                group=f"GeneTrek · {gs.category or 'Other'}",
                symbols=snap.matching_identifiers(members),
                n_genes=len(members),
                description=gs.description,
                secondary=not gs.is_primary and not used(gs.name),
            )
        )

    if data_dir is not None:
        for name, genes in sorted(load_genesets(data_dir).items()):
            sources.append(
                GenesetSource(
                    key=f"local:{name}",
                    label=name,
                    group=GROUP_LOCAL,
                    symbols=frozenset(genes),
                    n_genes=len(genes),
                )
            )
    return sources


def selected_symbols(
    sources: Iterable[GenesetSource], keys: Iterable[str]
) -> FrozenSet[str]:
    """The union of the selected sets."""
    wanted = set(keys)
    out: Set[str] = set()
    for source in sources:
        if source.key in wanted:
            out.update(source.symbols)
    return frozenset(out)


def geneset_filter_expr(column: str, symbols: Iterable[str]) -> pl.Expr:
    """True where the ``,``/``&`` joined gene cell names any of ``symbols``."""
    wanted = sorted(set(symbols))
    return (
        pl.col(column)
        .cast(pl.String)
        .str.to_uppercase()
        .str.extract_all(_GENE_TOKEN)
        .list.eval(pl.element().is_in(wanted))
        .list.any()
        .fill_null(False)
    )


def cell_matches(value: object, symbols: FrozenSet[str]) -> bool:
    """Python twin of :func:`geneset_filter_expr`, for row dicts."""
    if value is None:
        return False
    text = str(value).upper().replace("&", ",")
    return any(part.strip() in symbols for part in text.split(",") if part.strip())
