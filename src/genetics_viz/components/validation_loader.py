"""Shared validation data loading utilities."""

import csv
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from genetics_viz.utils.validation_badges import build_validation_badge

#: ``validations/snvs.tsv``.
SNV_VALIDATION_HEADER = [
    "FID", "Variant", "Sample", "User", "Inheritance", "Validation",
    "Comment", "Ignore", "Timestamp",
]  # fmt: skip

#: ``validations/svs.tsv``: the SNV columns plus the curated breakpoints.
SV_VALIDATION_HEADER = [
    "FID", "Variant", "Sample", "User", "Inheritance", "Validation",
    "CuratedStart", "CuratedEnd", "Comment", "Ignore", "Timestamp",
]  # fmt: skip


def validated_samples(
    validation_file: Path, family_id: str, variant_key: str
) -> Set[str]:
    """Samples of a family with any validation row for this variant."""
    if not validation_file.exists():
        return set()
    with open(validation_file, "r") as f:
        return {
            row.get("Sample") or ""
            for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
            if row.get("FID") == family_id and row.get("Variant") == variant_key
        }


def parents_to_validate(
    inheritance: str,
    parents: Dict[str, Optional[str]],
    validated: Set[str],
    both_parent_modes: Iterable[str] = ("either",),
) -> List[str]:
    """Parents to record as carriers when a child's call is validated present.

    ``maternal`` gives the mother, ``paternal`` the father, a mode in
    ``both_parent_modes`` both (mother first). A parent already validated for
    the variant, or not in the pedigree, is left out.
    """
    if inheritance == "maternal":
        wanted = [parents.get("mother")]
    elif inheritance == "paternal":
        wanted = [parents.get("father")]
    elif inheritance in both_parent_modes:
        wanted = [parents.get("mother"), parents.get("father")]
    else:
        wanted = []
    return [p for p in wanted if p and p not in validated]


def load_validation_map(
    validation_file_path: Path, family_id: str | None = None
) -> Dict[Tuple[str, str], List[Tuple[str, str, str, str, str, str, str]]]:
    """Load validation data from snvs.tsv or svs.tsv into a lookup map.

    Args:
        validation_file_path: Path to the validations file (snvs.tsv or svs.tsv)
        family_id: Optional family ID to filter by

    Returns:
        Dictionary mapping (variant_key, sample_id) to list of
        (validation_status, inheritance, comment, ignore, curated_start, curated_end, timestamp)
    """
    validation_map: Dict[
        Tuple[str, str], List[Tuple[str, str, str, str, str, str, str]]
    ] = {}

    if not validation_file_path.exists():
        return validation_map

    with open(validation_file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            fid = row.get("FID")
            variant_key = row.get("Variant")
            sample_id = row.get("Sample")
            validation_status = row.get("Validation")
            inheritance = row.get("Inheritance")
            comment = row.get("Comment", "")
            ignore = row.get("Ignore", "0")
            curated_start = row.get("CuratedStart", "")
            curated_end = row.get("CuratedEnd", "")
            timestamp = row.get("Timestamp", "")
            user = row.get("User", "")

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
                        curated_start or "",
                        curated_end or "",
                        timestamp or "",
                        user or "",
                    )
                )

    return validation_map


def add_validation_status_to_row(
    row: Dict[str, Any],
    validation_map: Dict[
        Tuple[str, str], List[Tuple[str, str, str, str, str, str, str]]
    ],
    variant_key: str,
    sample_id: str,
) -> None:
    """Add Validation and ValidationInheritance fields to a row.

    Args:
        row: The row dict to modify
        validation_map: Mapping from (variant_key, sample_id) to validations
        variant_key: The variant key (e.g. chr:pos:ref:alt or chr:start-end:type)
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
            row["Validation_badge"] = None
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
            # Check inheritance - prioritize de novo, then homozygous,
            # then first non-empty inheritance from present validations
            present = [v for v in validations if v[0] in ("present", "in phase MNV")]
            inh_values = [v[1] for v in present if v[1]]
            if "de novo" in inh_values:
                row["ValidationInheritance"] = "de novo"
            elif "homozygous" in inh_values:
                row["ValidationInheritance"] = "homozygous"
            elif inh_values:
                row["ValidationInheritance"] = inh_values[0]
            else:
                row["ValidationInheritance"] = ""
        elif "absent" in unique_validations:
            row["Validation"] = "absent"
            row["ValidationInheritance"] = ""
        else:
            row["Validation"] = "uncertain"
            row["ValidationInheritance"] = ""

        row["Validation_badge"] = build_validation_badge(
            row["Validation"], row["ValidationInheritance"], validations
        )
    else:
        row["Validation"] = ""
        row["ValidationInheritance"] = ""
        row["Validation_badge"] = None
