"""Pedigree consistency across the cohorts of one data directory.

The cohorts of a data directory share one ``families/`` tree, but each has its
own pedigree, and ghfc-ngs builds a family from the rows of the cohort it runs
(``FamilyPedigree.onDisk`` in ``main.nf``). Two cohorts listing one family ID
with other rows therefore overwrite each other's results. Barcodes have no
aliases, so one barcode under two family IDs is one family known by two codes
(cohort-specific IDs): it is analysed twice and its curation is split.

:func:`scan` reads every cohort pedigree as the pipeline does (positional
columns, the first row of each member, codes normalised by
``pedigree_edit``'s ``norm_*``) and returns :class:`Finding` s, worded in
``config/pedigree_consistency.yaml``. Curators acknowledge the ones that are
expected (``validations/pedigree_consistency.tsv``); a finding's key holds a
digest of the rows involved, so a change of rows raises it again. Nothing here
changes a pedigree.
"""

import csv
import hashlib
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, FrozenSet, Iterable, List, Mapping, Optional, Set, Tuple

import yaml

from ..models import is_dummy_id
from .auth import can_write
from .pedigree_edit import norm_parent, norm_phenotype, norm_sex
from .tsv import append_tsv_row

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "pedigree_consistency.yaml"
with open(_CONFIG_PATH, "r") as _f:
    CODES: Dict[str, Dict[str, str]] = dict((yaml.safe_load(_f) or {}).get("codes", {}))

SEVERITY_RANK = {"high": 3, "medium": 2, "info": 1}

#: ``(father, mother, sex, phenotype)``, normalised.
Row = Tuple[str, str, str, str]
#: ``(cohort, family ID)``.
Place = Tuple[str, str]


def code_info(code: str) -> Dict[str, str]:
    return CODES.get(
        code, {"severity": "info", "label": code, "short": code, "text": ""}
    )


# ---------------------------------------------------------------------------
# Reading the cohort pedigrees
# ---------------------------------------------------------------------------


def read_pedigree_rows(path: Path) -> Dict[str, Dict[str, Row]]:
    """``FID -> member -> Row``, as ghfc-ngs reads a cohort pedigree file."""
    with open(path, "r", encoding="utf-8", newline="") as f:
        return rows_from_lines(f)


def rows_from_lines(lines: Iterable[str]) -> Dict[str, Dict[str, Row]]:
    """``FID -> member -> Row`` from a pedigree's lines, as ghfc-ngs reads them.

    Positional columns (family, sample, father, mother, sex, phenotype), ``#``
    and blank lines skipped, a first line starting ``FID`` skipped, the first
    row of a member kept. A row of fewer than five columns is skipped.
    """
    out: Dict[str, Dict[str, Row]] = {}
    for index, line in enumerate(lines):
        if line.startswith("#") or not line.strip():
            continue
        cells = [c.strip() for c in line.rstrip("\r\n").split("\t")]
        if index == 0 and cells[0] == "FID":
            continue
        if len(cells) < 5 or not cells[0] or not cells[1]:
            continue
        members = out.setdefault(cells[0], {})
        if cells[1] in members:
            continue
        members[cells[1]] = (
            norm_parent(cells[2]),
            norm_parent(cells[3]),
            norm_sex(cells[4]),
            norm_phenotype(cells[5]) if len(cells) > 5 else "0",
        )
    return out


def _sequenced(members: Mapping[str, Row]) -> FrozenSet[str]:
    return frozenset(m for m in members if not is_dummy_id(m))


# ---------------------------------------------------------------------------
# Findings
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Finding:
    """One disagreement between pedigrees, about two or more families."""

    code: str
    #: Every ``(cohort, FID)`` concerned, sorted.
    places: Tuple[Place, ...]
    #: The barcodes concerned (shared ones for an alias).
    samples: Tuple[str, ...] = ()
    #: What differs, one line each (``S3: father 0 (PMS) / S1 (DEMO)``).
    details: Tuple[str, ...] = ()
    #: A digest of the rows involved: a change raises the finding again.
    digest: str = ""

    @property
    def info(self) -> Dict[str, str]:
        return code_info(self.code)

    @property
    def severity(self) -> str:
        return self.info.get("severity", "info")

    @property
    def rank(self) -> int:
        return SEVERITY_RANK.get(self.severity, 0)

    @property
    def family_ids(self) -> Tuple[str, ...]:
        return tuple(dict.fromkeys(fid for _, fid in self.places))

    @property
    def cohorts(self) -> Tuple[str, ...]:
        return tuple(dict.fromkeys(c for c, _ in self.places))

    @property
    def key(self) -> str:
        """What an acknowledgement is saved under."""
        where = ",".join(f"{c}:{f}" for c, f in self.places)
        return f"{self.code}|{where}|{self.digest}"

    @property
    def acknowledgeable(self) -> bool:
        return self.severity != "info"

    def others(self, family_id: str, cohort: Optional[str] = None) -> List[Place]:
        """The places other than ``(cohort, family_id)`` (any cohort when ``None``)."""
        return [
            (c, f)
            for c, f in self.places
            if not (f == family_id and (cohort is None or c == cohort))
        ]


def _digest(rows: Iterable[Tuple[Any, ...]]) -> str:
    text = "\n".join("\t".join(map(str, r)) for r in sorted(rows))
    return hashlib.sha1(text.encode()).hexdigest()[:10]


_FIELDS = ("father", "mother", "sex", "phenotype")


def _row_details(
    by_cohort: Mapping[str, Mapping[str, Row]],
    members: Iterable[str],
    fields: Tuple[int, ...],
) -> List[str]:
    """``S3: father 0 (PMS) / S1 (dataset_paper)``, per member and field that differ."""
    out: List[str] = []
    for member in sorted(members):
        parts = []
        for i in fields:
            values = {
                c: rows[member][i] for c, rows in by_cohort.items() if member in rows
            }
            if len(set(values.values())) > 1:
                parts.append(
                    f"{_FIELDS[i]} "
                    + " / ".join(f"{v} ({c})" for c, v in sorted(values.items()))
                )
        if parts:
            out.append(f"{member}: " + "; ".join(parts))
    return out


def _same_fid(fid: str, by_cohort: Dict[str, Dict[str, Row]]) -> Optional[Finding]:
    """The finding about one FID listed by several cohorts, or ``None``."""
    cohorts = sorted(by_cohort)
    places = tuple((c, fid) for c in cohorts)
    seq = {c: _sequenced(by_cohort[c]) for c in cohorts}
    digest = _digest((c, m, *r) for c in cohorts for m, r in by_cohort[c].items())
    if not frozenset.intersection(*seq.values()):
        return Finding(
            "FID_COLLISION",
            places,
            tuple(sorted(frozenset.union(*seq.values()))),
            tuple(f"{c}: {', '.join(sorted(seq[c]))}" for c in cohorts),
            digest,
        )
    if len(set(seq.values())) > 1:
        every = frozenset.union(*seq.values())
        details = [
            f"{c}: no {', '.join(sorted(every - seq[c]))}"
            for c in cohorts
            if every - seq[c]
        ] + _row_details(by_cohort, frozenset.intersection(*seq.values()), (0, 1, 2, 3))
        return Finding(
            "SAME_FID_MEMBERS_DIFFER",
            places,
            tuple(sorted(every)),
            tuple(details),
            digest,
        )
    members = set().union(*(by_cohort[c] for c in cohorts))
    # Dummies count here: the pipeline's rows hold them.
    shared = [m for m in members if all(m in by_cohort[c] for c in cohorts)]
    lacking = sorted(m for m in members if m not in shared)
    rows = _row_details(by_cohort, shared, (0, 1, 2))
    if rows or lacking:
        rows += [
            f"{m}: only in {', '.join(c for c in cohorts if m in by_cohort[c])}"
            for m in lacking
        ]
        return Finding(
            "SAME_FID_ROWS_DIFFER",
            places,
            tuple(sorted(seq[cohorts[0]])),
            tuple(rows),
            digest,
        )
    pheno = _row_details(by_cohort, shared, (3,))
    if pheno:
        return Finding(
            "SAME_FID_PHENOTYPE_DIFFERS",
            places,
            tuple(sorted(seq[cohorts[0]])),
            tuple(pheno),
            digest,
        )
    return None


def _aliases(
    pedigrees: Mapping[str, Mapping[str, Mapping[str, Row]]],
) -> List[Finding]:
    """Families with different IDs sharing a sequenced barcode."""
    where: Dict[str, Set[Place]] = {}
    for cohort, families in pedigrees.items():
        for fid, members in families.items():
            for sample in _sequenced(members):
                where.setdefault(sample, set()).add((cohort, fid))
    # Each pair of family IDs sharing a barcode, with the places of each.
    pairs: Dict[Tuple[str, str], Set[str]] = {}
    for sample, places in where.items():
        fids = sorted({f for _, f in places})
        for i, a in enumerate(fids):
            for b in fids[i + 1 :]:
                pairs.setdefault((a, b), set()).add(sample)
    out: List[Finding] = []
    for (a, b), shared in sorted(pairs.items()):
        places_a = sorted(p for s in shared for p in where[s] if p[1] == a)
        places_b = sorted(p for s in shared for p in where[s] if p[1] == b)
        places = tuple(sorted(set(places_a) | set(places_b)))

        def members(places_: Iterable[Place]) -> Dict[str, Row]:
            merged: Dict[str, Row] = {}
            for c, f in places_:
                for m, r in pedigrees[c][f].items():
                    merged.setdefault(m, r)
            return merged

        rows_a, rows_b = members(set(places_a)), members(set(places_b))
        seq_a, seq_b = _sequenced(rows_a), _sequenced(rows_b)
        if seq_a == seq_b:
            code = "ALIAS_SAME_MEMBERS"
        elif seq_a <= seq_b or seq_b <= seq_a:
            code = "ALIAS_SUBSET"
        else:
            code = "ALIAS_OVERLAP"
        details = [
            f"{a}: {len(seq_a)} sequenced, {b}: {len(seq_b)}; shared: "
            + ", ".join(sorted(shared))
        ]
        # A shared sample with another sex, or another parent that both
        # families hold, contradicts itself: that outranks the alias.
        conflicts: List[str] = []
        for s in sorted(shared):
            ra, rb = rows_a[s], rows_b[s]
            if ra[2] != rb[2] and "0" not in (ra[2], rb[2]):
                conflicts.append(f"{s}: sex {ra[2]} ({a}) / {rb[2]} ({b})")
            # Only parents both families hold can be compared (or none): a
            # parent one family lacks is the subset finding's matter.
            known = (seq_a & seq_b) | {"0"}
            for i, role in ((0, "father"), (1, "mother")):
                pa, pb = ra[i], rb[i]
                if pa != pb and pa in known and pb in known:
                    conflicts.append(f"{s}: {role} {pa} ({a}) / {pb} ({b})")
        digest = _digest(
            [
                (f, m, *r)
                for f, rows in ((a, rows_a), (b, rows_b))
                for m, r in rows.items()
            ]
        )
        if conflicts:
            out.append(
                Finding(
                    "BARCODE_ROWS_DIFFER",
                    places,
                    tuple(sorted(shared)),
                    tuple(conflicts + details),
                    digest,
                )
            )
        out.append(Finding(code, places, tuple(sorted(shared)), tuple(details), digest))
    return out


def scan_pedigrees(
    pedigrees: Mapping[str, Mapping[str, Mapping[str, Row]]],
) -> List[Finding]:
    """Every finding over ``cohort -> FID -> member -> Row``, worst first."""
    by_fid: Dict[str, Dict[str, Dict[str, Row]]] = {}
    for cohort, families in pedigrees.items():
        for fid, members in families.items():
            by_fid.setdefault(fid, {})[cohort] = dict(members)
    out: List[Finding] = []
    for fid, by_cohort in sorted(by_fid.items()):
        if len(by_cohort) > 1:
            finding = _same_fid(fid, by_cohort)
            if finding is not None:
                out.append(finding)
    out.extend(_aliases(pedigrees))
    return sorted(out, key=lambda f: (-f.rank, f.code, f.places))


@dataclass
class Consistency:
    """A data directory's findings, looked up by family."""

    findings: List[Finding] = field(default_factory=list)

    def for_family(self, family_id: str) -> List[Finding]:
        return [f for f in self.findings if family_id in f.family_ids]

    def alias_groups(self) -> Dict[str, FrozenSet[str]]:
        """``FID -> every FID known to be the same family`` (shared barcodes).

        Only same-members and subset aliases join: an overlap may be a
        misplaced member. A family with no alias is absent.
        """
        parent: Dict[str, str] = {}

        def find(x: str) -> str:
            while parent.setdefault(x, x) != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for f in self.findings:
            if f.code in ("ALIAS_SAME_MEMBERS", "ALIAS_SUBSET"):
                ids = f.family_ids
                for other in ids[1:]:
                    parent[find(other)] = find(ids[0])
        groups: Dict[str, Set[str]] = {}
        for fid in parent:
            groups.setdefault(find(fid), set()).add(fid)
        return {fid: frozenset(g) for g in groups.values() if len(g) > 1 for fid in g}


#: ``(path, mtime, size)`` of every cohort pedigree -> its consistency.
_cache: Dict[Tuple[Tuple[str, float, int], ...], Consistency] = {}


def scan(store: Any) -> Consistency:
    """Blocking: the store's cohort pedigrees, read and compared.

    Cached on the pedigree files' stats. Never raises: a file that cannot be
    read is left out and logged.
    """
    files: List[Tuple[str, Path]] = []
    for name, cohort in sorted(store.cohorts.items()):
        path = getattr(cohort, "pedigree_file", None)
        if path is not None:
            files.append((name, Path(path)))
    stats = []
    for _name, path in files:
        try:
            st = path.stat()
            stats.append((str(path), st.st_mtime, st.st_size))
        except OSError:
            continue
    key = tuple(stats)
    if key in _cache:
        return _cache[key]
    pedigrees: Dict[str, Dict[str, Dict[str, Row]]] = {}
    for name, path in files:
        try:
            pedigrees[name] = read_pedigree_rows(path)
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
    result = Consistency(scan_pedigrees(pedigrees))
    _cache.clear()
    _cache[key] = result
    return result


def clear_consistency_cache() -> None:
    _cache.clear()
    _ack_cache.clear()


# ---------------------------------------------------------------------------
# Acknowledgements
# ---------------------------------------------------------------------------

ACK_HEADER = ["Key", "Code", "Families", "Status", "Note", "User", "Timestamp"]
EXPLAINED = "explained"
WITHDRAWN = "withdrawn"


@dataclass(frozen=True)
class Ack:
    status: str
    note: str
    user: str
    timestamp: str

    @property
    def date(self) -> str:
        return self.timestamp[:10]


def ack_file(data_dir: Path) -> Path:
    return Path(data_dir) / "validations" / "pedigree_consistency.tsv"


def save_ack(
    data_dir: Path,
    finding: Finding,
    note: str,
    *,
    user: str,
    status: str = EXPLAINED,
) -> None:
    """Record that ``finding`` is expected (``withdrawn`` takes it back)."""
    if not can_write():
        raise PermissionError("Only curators and administrators can explain a finding")
    if not user:
        raise ValueError("No user to record the explanation under")
    if status not in (EXPLAINED, WITHDRAWN):
        raise ValueError(f"Unknown status: {status}")
    append_tsv_row(
        ack_file(data_dir),
        ACK_HEADER,
        [
            finding.key,
            finding.code,
            ",".join(f"{c}:{f}" for c, f in finding.places),
            status,
            " ".join((note or "").split()),
            user,
            datetime.now().isoformat(timespec="seconds"),
        ],
    )


_ack_cache: Dict[Tuple[str, float, int], Dict[str, Ack]] = {}


def load_acks(data_dir: Path) -> Dict[str, Ack]:
    """``finding key -> standing explanation`` (withdrawn ones left out). Never raises."""
    path = ack_file(data_dir)
    try:
        st = path.stat()
    except OSError:
        return {}
    key = (str(path), st.st_mtime, st.st_size)
    if key not in _ack_cache:
        latest: Dict[str, Ack] = {}
        try:
            with open(path, "r", newline="") as f:
                for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
                    if row.get("Key"):
                        latest[row["Key"]] = Ack(
                            row.get("Status") or "",
                            row.get("Note") or "",
                            row.get("User") or "",
                            row.get("Timestamp") or "",
                        )
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            return {}
        _ack_cache.clear()
        _ack_cache[key] = {k: a for k, a in latest.items() if a.status == EXPLAINED}
    return _ack_cache[key]


def open_findings(
    findings: Iterable[Finding], acks: Mapping[str, Ack]
) -> List[Finding]:
    """The findings to review: not info, not explained."""
    return [f for f in findings if f.acknowledgeable and f.key not in acks]
