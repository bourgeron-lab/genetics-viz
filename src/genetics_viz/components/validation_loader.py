"""Shared validation data loading utilities."""

import csv
from pathlib import Path
from typing import Any, Dict, List, Tuple


def load_validation_map(
    validation_file_path: Path, family_id: str | None = None
) -> Dict[Tuple[str, str], List[Tuple[str, str]]]:
    """Load validation data from snvs.tsv into a lookup map.

    Args:
        validation_file_path: Path to the validations/snvs.tsv file
        family_id: Optional family ID to filter by

    Returns:
        Dictionary mapping (variant_key, sample_id) to list of (validation_status, inheritance)
    """
    validation_map: Dict[Tuple[str, str], List[Tuple[str, str]]] = {}

    if not validation_file_path.exists():
        return validation_map

    with open(validation_file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fid = row.get("FID")
            variant_key = row.get("Variant")
            sample_id = row.get("Sample")
            validation_status = row.get("Validation")
            inheritance = row.get("Inheritance")

            # Filter by family_id if provided
            if family_id is not None and fid != family_id:
                continue

            if variant_key and sample_id:
                map_key = (variant_key, sample_id)
                if map_key not in validation_map:
                    validation_map[map_key] = []
                validation_map[map_key].append(
                    (validation_status or "", inheritance or "")
                )

    return validation_map


def add_validation_status_to_row(
    row: Dict[str, Any],
    validation_map: Dict[Tuple[str, str], List[Tuple[str, str]]],
    variant_key: str,
    sample_id: str,
) -> None:
    """Add Validation and ValidationInheritance fields to a row.

    Args:
        row: The row dict to modify
        validation_map: Mapping from (variant_key, sample_id) to validations
        variant_key: The variant key (e.g. chr:pos:ref:alt)
        sample_id: The sample ID
    """
    map_key = (variant_key, sample_id)

    if map_key in validation_map:
        validations = validation_map[map_key]
        validation_statuses = [v[0] for v in validations]
        unique_validations = set(validation_statuses)

        if len(unique_validations) > 1:
            row["Validation"] = "conflicting"
            row["ValidationInheritance"] = ""
        elif "present" in unique_validations:
            row["Validation"] = "present"
            is_de_novo = any(
                v[1] == "de novo" for v in validations if v[0] == "present"
            )
            row["ValidationInheritance"] = "de novo" if is_de_novo else ""
        elif "absent" in unique_validations:
            row["Validation"] = "absent"
            row["ValidationInheritance"] = ""
        else:
            row["Validation"] = "uncertain"
            row["ValidationInheritance"] = ""
    else:
        row["Validation"] = ""
        row["ValidationInheritance"] = ""
