"""Geneset loading utilities.

Loads gene sets from ``<data_dir>/params/genesets/*.tsv`` files.
Each TSV has one gene symbol per line (first line is a header).
"""

from pathlib import Path
from typing import Dict, Set


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
