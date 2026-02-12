"""Gene scoring and color coding based on genesets."""

import csv
from pathlib import Path
from typing import Dict, List, Tuple

import yaml


class GeneScorer:
    """Singleton class for gene scoring based on genesets."""

    _instance = None
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if not self._initialized:
            self._load_data()
            GeneScorer._initialized = True

    def _load_data(self):
        """Load genesets configuration and data."""
        config_dir = Path(__file__).parent.parent / "config"

        # Load genesets.yaml for scores
        yaml_path = config_dir / "genesets.yaml"
        with open(yaml_path, "r") as f:
            self.genesets_config = yaml.safe_load(f)

        # Load genesets.tsv
        tsv_path = config_dir / "genesets.tsv"
        self.gene_data: Dict[str, Dict] = {}  # symbol -> gene info
        self.ensg_to_symbol: Dict[str, str] = {}  # ENSG -> symbol

        with open(tsv_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                symbol = row.get("symbol", "")
                ensg = row.get("ensembl_gene_id", "")

                if symbol:
                    self.gene_data[symbol] = row
                    if ensg:
                        self.ensg_to_symbol[ensg] = symbol

    def get_gene_score_and_sets(self, gene_identifier: str) -> Tuple[float, List[str]]:
        """Get score and geneset list for a gene symbol or ENSG ID.

        Args:
            gene_identifier: Gene symbol (e.g., 'BRCA1') or ENSG ID (e.g., 'ENSG00000012048')

        Returns:
            Tuple of (score, list of geneset names)
        """
        # Check if it's an ENSG ID
        if gene_identifier.startswith("ENSG"):
            gene_identifier = self.ensg_to_symbol.get(gene_identifier, gene_identifier)

        gene_row = self.gene_data.get(gene_identifier)
        if not gene_row:
            return (0.0, [])

        score = 0.0
        geneset_names = []

        for geneset_name, geneset_config in self.genesets_config.items():
            # Check if gene is in this geneset
            if gene_row.get(geneset_name, "False") == "True":
                score += geneset_config.get("score", 0)
                geneset_names.append(geneset_name)

        return (score, geneset_names)

    def get_gene_color(self, gene_identifier: str) -> str:
        """Get color for a gene based on its score.

        Args:
            gene_identifier: Gene symbol or ENSG ID

        Returns:
            Hex color string (white to dark red gradient)
        """
        score, _ = self.get_gene_score_and_sets(gene_identifier)
        return self.score_to_color(score)

    @staticmethod
    def score_to_color(score: float) -> str:
        """Convert score to color (linear gradient white to dark red).

        Args:
            score: Gene score (0 = white, 8+ = dark red)

        Returns:
            Hex color string
        """
        # Clamp score between 0 and 8
        clamped_score = max(0, min(8, score))

        # Linear interpolation from white (255,255,255) to dark red (139,0,0)
        ratio = clamped_score / 8.0

        # RGB values
        r = int(255 - (255 - 139) * ratio)
        g = int(255 * (1 - ratio))
        b = int(255 * (1 - ratio))

        return f"#{r:02x}{g:02x}{b:02x}"

    def get_gene_tooltip(self, gene_identifier: str) -> str:
        """Get tooltip text for a gene showing its genesets.

        Args:
            gene_identifier: Gene symbol or ENSG ID

        Returns:
            Tooltip string with geneset names
        """
        score, geneset_names = self.get_gene_score_and_sets(gene_identifier)

        if not geneset_names:
            return gene_identifier

        return f"{gene_identifier}: {', '.join(geneset_names)}"

    def reload(self):
        """Reload gene scoring data from config files."""
        self._load_data()


# Global instance
_gene_scorer = None


def get_gene_scorer() -> GeneScorer:
    """Get the global GeneScorer instance."""
    global _gene_scorer
    if _gene_scorer is None:
        _gene_scorer = GeneScorer()
    return _gene_scorer


def reload_gene_scoring():
    """Reload gene scoring configuration."""
    scorer = get_gene_scorer()
    scorer.reload()
