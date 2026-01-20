"""Shared validation data loading utilities."""

import csv
from pathlib import Path
from typing import Any, Dict, List, Tuple


def load_validation_map(
    validation_file_path: Path, family_id: str | None = None
) -> Dict[Tuple[str, str], List[Tuple[str, str, str, str]]]:
    """Load validation data from snvs.tsv into a lookup map.

    Args:
        validation_file_path: Path to the validations/snvs.tsv file
        family_id: Optional family ID to filter by

    Returns:
        Dictionary mapping (variant_key, sample_id) to list of
        (validation_status, inheritance, comment, ignore)
    """
    validation_map: Dict[Tuple[str, str], List[Tuple[str, str, str, str]]] = {}

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
            comment = row.get("Comment", "")
            ignore = row.get("Ignore", "0")

            # Filter by family_id if provided
            if family_id is not None and fid != family_id:
                continue

            if variant_key and sample_id:
                map_key = (variant_key, sample_id)
                if map_key not in validation_map:
                    validation_map[map_key] = []
                validation_map[map_key].append(
                    (
                        validation_status or "",
                        inheritance or "",
                        comment or "",
                        ignore or "0",
                    )
                )

    return validation_map


def add_validation_status_to_row(
    row: Dict[str, Any],
    validation_map: Dict[Tuple[str, str], List[Tuple[str, str, str, str]]],
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
        all_validations = validation_map[map_key]
        # Filter out ignored validations (Ignore=1)
        validations = [v for v in all_validations if v[3] != "1"]

        if not validations:
            # All validations are ignored
            row["Validation"] = ""
            row["ValidationInheritance"] = ""
            return

        validation_statuses = [v[0] for v in validations]
        # Normalize "in phase MNV" to "present" for conflict detection
        normalized_statuses = [
            "present" if s == "in phase MNV" else s for s in validation_statuses
        ]
        unique_validations = set(normalized_statuses)

        if len(unique_validations) > 1:
            row["Validation"] = "conflicting"
            row["ValidationInheritance"] = ""
        elif "present" in unique_validations:
            # Check if any is specifically "in phase MNV"
            if "in phase MNV" in validation_statuses:
                row["Validation"] = "in phase MNV"
            else:
                row["Validation"] = "present"
            # Check inheritance - prioritize de novo, then homozygous
            is_de_novo = any(
                v[1] == "de novo"
                for v in validations
                if v[0] in ("present", "in phase MNV")
            )
            is_homozygous = any(
                v[1] == "homozygous"
                for v in validations
                if v[0] in ("present", "in phase MNV")
            )
            if is_de_novo:
                row["ValidationInheritance"] = "de novo"
            elif is_homozygous:
                row["ValidationInheritance"] = "homozygous"
            else:
                row["ValidationInheritance"] = ""
        elif "absent" in unique_validations:
            row["Validation"] = "absent"
            row["ValidationInheritance"] = ""
        else:
            row["Validation"] = "uncertain"
            row["ValidationInheritance"] = ""
    else:
        row["Validation"] = ""
        row["ValidationInheritance"] = ""
