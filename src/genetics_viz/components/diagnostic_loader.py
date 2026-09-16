"""Shared diagnostic data loading utilities."""

import csv
import fcntl
from pathlib import Path
from typing import Any, Dict, List, Tuple

from genetics_viz.utils.diagnostic_badges import build_diagnostic_badge

# Verdict used when curators disagree on the same variant for the same sample.
CONFLICTING_DIAGNOSTIC = "conflicting"

# TSV header for diagnostics files
DIAGNOSTIC_HEADER = (
    "FID\tVariant\tGene\tImpact\tSample\tUser\tTimestamp\tComment\tIgnore\tDiagnostic\n"
)


def ensure_diagnostic_file(file_path: Path) -> None:
    """Create the diagnostics file with header if it doesn't exist."""
    file_path.parent.mkdir(parents=True, exist_ok=True)
    if not file_path.exists():
        with open(file_path, "w") as f:
            f.write(DIAGNOSTIC_HEADER)


def load_diagnostic_map(
    diagnostic_file_path: Path, family_id: str | None = None
) -> Dict[Tuple[str, str], List[Tuple[str, str, str, str, str]]]:
    """Load diagnostic data from snvs.tsv or svs.tsv into a lookup map.

    Args:
        diagnostic_file_path: Path to the diagnostics file
        family_id: Optional family ID to filter by

    Returns:
        Dictionary mapping (variant_key, sample_id) to list of
        (diagnostic_value, user, timestamp, comment, ignore)
    """
    diagnostic_map: Dict[Tuple[str, str], List[Tuple[str, str, str, str, str]]] = {}

    if not diagnostic_file_path.exists():
        return diagnostic_map

    with open(diagnostic_file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fid = row.get("FID")
            variant_key = row.get("Variant")
            sample_id = row.get("Sample")
            diagnostic = row.get("Diagnostic", "")
            user = row.get("User", "")
            timestamp = row.get("Timestamp", "")
            comment = row.get("Comment", "")
            ignore = row.get("Ignore", "0")

            if family_id is not None and fid != family_id:
                continue

            if variant_key and sample_id:
                map_key = (variant_key, sample_id)
                if map_key not in diagnostic_map:
                    diagnostic_map[map_key] = []
                diagnostic_map[map_key].append(
                    (diagnostic, user, timestamp, comment, ignore)
                )

    return diagnostic_map


def add_diagnostic_status_to_row(
    row: Dict[str, Any],
    diagnostic_map: Dict[Tuple[str, str], List[Tuple[str, str, str, str, str]]],
    variant_key: str,
    sample_id: str,
) -> None:
    """Add Diagnostic and Diagnostic_badge fields to a row.

    Args:
        row: The row dict to modify
        diagnostic_map: Mapping from (variant_key, sample_id) to diagnostics
        variant_key: The variant key
        sample_id: The sample ID
    """
    map_key = (variant_key, sample_id)

    if map_key in diagnostic_map:
        all_diagnostics = diagnostic_map[map_key]
        # Filter out ignored entries
        diagnostics = [d for d in all_diagnostics if d[4] != "1"]

        if not diagnostics:
            row["Diagnostic"] = ""
            row["Diagnostic_badge"] = None
            return

        diagnostic_values = [d[0] for d in diagnostics]
        unique_diagnostics = set(diagnostic_values)

        if len(unique_diagnostics) > 1:
            row["Diagnostic"] = CONFLICTING_DIAGNOSTIC
        else:
            row["Diagnostic"] = diagnostic_values[0]

        # Build badge with non-ignored entries for tooltip
        badge_data = [(d[0], d[1], d[2], d[3]) for d in diagnostics]
        row["Diagnostic_badge"] = build_diagnostic_badge(row["Diagnostic"], badge_data)
    else:
        row["Diagnostic"] = ""
        row["Diagnostic_badge"] = None


def save_diagnostic_entry(
    diagnostic_file: Path,
    family_id: str,
    variant_key: str,
    gene: str,
    impact: str,
    sample: str,
    user: str,
    timestamp: str,
    comment: str,
    diagnostic: str,
) -> None:
    """Append a diagnostic entry to the TSV file with file locking."""
    ensure_diagnostic_file(diagnostic_file)
    with open(diagnostic_file, "a") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            f.write(
                f"{family_id}\t{variant_key}\t{gene}\t{impact}\t{sample}\t"
                f"{user}\t{timestamp}\t{comment}\t0\t{diagnostic}\n"
            )
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def update_diagnostic_ignore_status(
    diagnostic_file: Path,
    family_id: str,
    variant_key: str,
    sample: str,
    timestamp: str,
    ignore_value: str,
) -> bool:
    """Update the Ignore status for a specific diagnostic row.

    Uses read-all, modify, write-back pattern with file locking.

    Returns:
        True if update was successful, False otherwise
    """
    if not diagnostic_file.exists():
        return False

    rows = []
    fieldnames: list[str] = []
    with open(diagnostic_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            rows.append(row)

    updated = False
    for row in rows:
        if (
            row.get("FID") == family_id
            and row.get("Variant") == variant_key
            and row.get("Sample") == sample
            and row.get("Timestamp") == timestamp
        ):
            row["Ignore"] = ignore_value
            updated = True
            break

    if not updated:
        return False

    with open(diagnostic_file, "w") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)

    return True


def load_family_diagnostics(
    snv_file: Path,
    sv_file: Path,
    family_id: str,
    sample_ids: List[str],
) -> List[Dict[str, Any]]:
    """Load the diagnostic entries for a family and set of samples.

    Rows are merged twice over, so a variant is listed once per verdict rather
    than once per saved record:

    * Curators first. Several curators recording the same variant for the same
      sample are one assessment - ``"conflicting"`` when they disagree, which is
      the aggregation :func:`add_diagnostic_status_to_row` already applies to
      the variant tables.
    * Then individuals. Samples whose assessment reached the *same* verdict for
      a variant fold into one row. Samples that reached different verdicts stay
      on separate rows, because a variant can legitimately be pathogenic in an
      affected child and benign in a parent, and collapsing that would report a
      per-individual assessment as a curator disagreement.

    ``Sample`` therefore lists every individual the row covers, in the pedigree
    order of ``sample_ids``, and ``User`` every curator involved, newest first.
    ``_entries`` keeps the underlying records so a caller can show who recorded
    what - see :func:`format_diagnostic_contributors`.

    Returns one row dict per (source, variant, verdict).
    """
    per_sample: Dict[Tuple[str, str, str], List[Dict[str, str]]] = {}
    sample_set = set(sample_ids)

    for diag_file in [snv_file, sv_file]:
        if not diag_file.exists():
            continue
        source = "snv" if diag_file == snv_file else "sv"
        with open(diag_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if (
                    row.get("FID") == family_id
                    and row.get("Sample") in sample_set
                    and row.get("Ignore", "0") != "1"
                ):
                    key = (source, row.get("Variant", ""), row.get("Sample", ""))
                    per_sample.setdefault(key, []).append(row)

    # Settle each sample's verdict, then key by it so same-verdict samples meet.
    merged: Dict[Tuple[str, str, str], List[Dict[str, str]]] = {}
    for (source, variant, _sample), rows in per_sample.items():
        verdicts = {r.get("Diagnostic", "") for r in rows}
        verdict = next(iter(verdicts)) if len(verdicts) == 1 else CONFLICTING_DIAGNOSTIC
        merged.setdefault((source, variant, verdict), []).extend(rows)

    sample_order = {sid: i for i, sid in enumerate(sample_ids)}
    entries: List[Dict[str, Any]] = []
    for (source, _variant, verdict), rows in merged.items():
        rows = sorted(rows, key=lambda r: r.get("Timestamp", ""), reverse=True)
        # dict.fromkeys de-duplicates while preserving order, for the curator
        # who recorded one variant twice and the sample carrying two records.
        users = list(dict.fromkeys(r.get("User", "") for r in rows if r.get("User")))
        samples = sorted(
            dict.fromkeys(r.get("Sample", "") for r in rows if r.get("Sample")),
            key=lambda s: sample_order.get(s, len(sample_order)),
        )
        entries.append(
            {
                **rows[0],
                "_source": source,
                "_entries": rows,
                "Diagnostic": verdict,
                "User": ", ".join(users),
                "Sample": ", ".join(samples),
            }
        )

    return entries


def format_diagnostic_contributors(entry: Dict[str, Any]) -> str:
    """Summarise who recorded what for a merged :func:`load_family_diagnostics` row.

    Reads the ``_entries`` the row was merged from and returns a line such as
    ``"alice: pathogenic (2026-03-01) - bob: benign (2026-04-02)"``, which tells
    the reader who disagreed when the row is marked ``"conflicting"``. A row
    covering more than one individual names the sample each record was made
    against, since the curator alone no longer identifies it.
    """
    records = entry.get("_entries", [])
    show_sample = len({r.get("Sample", "") for r in records}) > 1
    parts = []
    for row in records:
        who = row.get("User", "") or "unknown"
        if show_sample:
            who = f"{who} on {row.get('Sample', '?')}"
        verdict = row.get("Diagnostic", "") or "?"
        date = row.get("Timestamp", "").split(" ")[0].split("T")[0]
        parts.append(f"{who}: {verdict} ({date})" if date else f"{who}: {verdict}")
    return " - ".join(parts)
