"""Sample QC exports, shared by the family and cohort tabs: flags, samples, pairs.

Three tables, each worded for the user's role as the tabs are:

- ``flags``: one row per flag, with the curators' acknowledgement in its own
  columns (status, note, user, date). Acknowledged flags are always listed.
- ``samples``: sex, karyotype and quality metrics per sample, with the flags
  still counting and, apart, the acknowledged ones.
- ``pairs``: every pair tested, against the pedigree, with the pair's flags.

Readers get ``reader_label``, no evidence, an "Atypical" karyotype, no
observed relationship, only the pedigree-declared pairs, and an
acknowledgement's note only for the categories marked ``reader_note``.
"""

from typing import Any, Dict, Iterable, List, Optional, Tuple

from .export import PERSON_COLUMNS
from .sample_qc import (
    FlagRow,
    flag_title,
    flags_by_key,
    is_typical_karyotype,
    severity_info,
)
from .sample_qc_acks import Ack, AckKey, note_visible, split_acknowledged

Key = Tuple[str, str]
Table = Tuple[List[Dict[str, Any]], List[str]]


def _flag_texts(
    flags: Iterable[FlagRow], acks: Dict[AckKey, Ack], show_detail: bool
) -> Tuple[str, str]:
    """``(still counting, acknowledged)``: titles ``"; "``-joined, each once."""
    active, acked = split_acknowledged(flags, acks)
    counting = dict.fromkeys(flag_title(f, show_detail) for f in active)
    reviewed = dict.fromkeys(
        f"{flag_title(f, show_detail)} ({acks[f.ack_key].label.lower()})" for f in acked
    )
    return "; ".join(counting), "; ".join(reviewed)


def flags_export(
    flags: Iterable[FlagRow], acks: Dict[AckKey, Ack], show_detail: bool
) -> Table:
    """One row per flag, most severe first, acknowledged ones after."""
    columns = ["level", "FID", "sample_a", "sample_b", "code", "issue", "category",
               "severity"]  # fmt: skip
    if show_detail:
        columns.append("evidence")
    columns += ["ack_status", "ack_note", "ack_user", "ack_date"]
    active, acked = split_acknowledged(flags, acks)
    rows = []
    for f in active + acked:
        ack = acks.get(f.ack_key)
        rows.append(
            {
                "level": f.level,
                "FID": f.family_id,
                "sample_a": f.sample_a,
                "sample_b": f.sample_b or "",
                "code": f.code,
                "issue": flag_title(f, show_detail),
                "category": f.info.category_label,
                "severity": severity_info(f.severity)["key"],
                "evidence": f.evidence,
                "ack_status": ack.status if ack else "",
                "ack_note": (
                    (ack.note if note_visible(f, show_detail) else "(hidden)")
                    if ack and ack.note
                    else ""
                ),
                "ack_user": ack.user if ack else "",
                "ack_date": ack.date if ack else "",
            }
        )
    return rows, columns


#: ``samples.tsv`` columns written as they are, after the sex calls.
_SAMPLE_METRICS = [
    "x_copies", "y_copies", "somalier_x_copies", "somalier_y_copies",
    "x_het_rate", "x_sites", "depth_mean", "n_sites_called", "het_hom_alt_ratio",
    "contamination_charr", "charr_sites", "somalier_charr", "p_middling_ab",
    "trio_informative_sites", "trio_mendelian_error_rate",
]  # fmt: skip


def samples_export(
    sample_rows: Iterable[Dict[str, Any]],
    people: Dict[Key, Dict[str, Any]],
    flags: Iterable[FlagRow],
    acks: Dict[AckKey, Ack],
    show_detail: bool,
) -> Table:
    """One row per checked person: sex, karyotype, quality, then their flags.

    ``people`` maps ``(family_id, sample)`` to an ``export.person`` dict; rows
    of anyone else are left out.
    """
    columns = [*PERSON_COLUMNS, "pedigree_sex", "inferred_sex", "karyotype",
               "sex_evidence", *_SAMPLE_METRICS, "worst_severity", "flags",
               "acknowledged_flags"]  # fmt: skip
    by_key = flags_by_key(flags)
    rows = []
    for r in sample_rows:
        key = (r.get("family_id") or "", r.get("sample") or "")
        who = people.get(key)
        if who is None:
            continue
        karyotype = r.get("karyotype") or "unknown"
        if not show_detail and not is_typical_karyotype(karyotype):
            karyotype = "Atypical" if karyotype != "unknown" else "unknown"
        mine = by_key.get(key, [])
        active, _ = split_acknowledged(mine, acks)
        counting, reviewed = _flag_texts(mine, acks, show_detail)
        rows.append(
            {
                **who,
                "pedigree_sex": r.get("pedigree_sex"),
                "inferred_sex": r.get("inferred_sex"),
                "karyotype": karyotype,
                "sex_evidence": r.get("sex_evidence"),
                **{m: r.get(m) for m in _SAMPLE_METRICS},
                "worst_severity": (
                    severity_info(max(active, key=lambda f: f.rank).severity)["key"]
                    if active
                    else ""
                ),
                "flags": counting,
                "acknowledged_flags": reviewed,
            }
        )
    return rows, columns


#: ``pairs.tsv`` columns written as they are, after the relationships.
_PAIR_METRICS = [
    "king_kinship", "expected_kinship", "somalier_relatedness", "ibs0", "ibs2",
    "ibs0_ibs2_ratio", "hom_concordance", "n_sites", "x_ibs0", "x_ibs2",
    "po_informative_sites", "po_opposite_hom_rate",
]  # fmt: skip


def _declared(row: Dict[str, Any]) -> bool:
    same = (row.get("family_id_b") or row.get("family_id_a")) == row.get("family_id_a")
    return same and (row.get("expected_relationship") or "unrelated") != "unrelated"


def pairs_export(
    pair_rows: Iterable[Dict[str, Any]],
    flags: Iterable[FlagRow],
    acks: Dict[AckKey, Ack],
    show_detail: bool,
    first_family: Optional[str] = None,
) -> Table:
    """One row per pair (``pairs.tsv`` dicts), with the flags raised on it.

    Readers get the pedigree-declared pairs only. ``first_family`` puts that
    family's member first in a cross-family pair.
    """
    columns = ["FID_a", "sample_a", "FID_b", "sample_b", "expected_relationship"]
    if show_detail:
        columns.append("observed_relationship")
    columns += ["status", *_PAIR_METRICS, "flags", "acknowledged_flags"]
    by_pair: Dict[frozenset, List[FlagRow]] = {}
    for f in flags:
        if f.sample_b:
            by_pair.setdefault(frozenset(f.keys), []).append(f)
    rows = []
    for r in pair_rows:
        if not r.get("sample_a") or not r.get("sample_b"):
            continue
        if not show_detail and not _declared(r):
            continue
        a = (r.get("family_id_a") or "", r["sample_a"])
        b = (r.get("family_id_b") or a[0], r["sample_b"])
        if first_family is not None and b[0] == first_family and a[0] != first_family:
            a, b = b, a
        counting, reviewed = _flag_texts(
            by_pair.get(frozenset((a, b)), []), acks, show_detail
        )
        rows.append(
            {
                "FID_a": a[0],
                "sample_a": a[1],
                "FID_b": b[0],
                "sample_b": b[1],
                "expected_relationship": r.get("expected_relationship"),
                "observed_relationship": r.get("observed_relationship"),
                "status": r.get("status"),
                **{m: r.get(m) for m in _PAIR_METRICS},
                "flags": counting,
                "acknowledged_flags": reviewed,
            }
        )
    return rows, columns
