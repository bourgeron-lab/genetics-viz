"""Login backoff: per client address and per username, in memory.

After a few free failures, each further failure doubles the wait before the
next attempt, up to ``MAX_DELAY``. An attempt during the wait is refused
without checking the password, so it tells the client nothing. Unknown
usernames are counted like real ones, so a refusal does not reveal which
accounts exist.

The address is the one the server sees. Behind a reverse proxy that is the
proxy's, unless uvicorn trusts its forwarded headers (``FORWARDED_ALLOW_IPS``);
the per-address limit is set high for that reason, and the per-user one still
holds.

State is per process and lost on restart, which is enough to make online
guessing slow.
"""

import time
from dataclasses import dataclass
from typing import Callable

#: Failures allowed before the backoff starts, per username and per address.
FREE_FAILURES_USER = 5
FREE_FAILURES_ADDRESS = 20

BASE_DELAY = 2.0  # seconds, for the first failure past the free ones
MAX_DELAY = 15 * 60.0

#: A key with no failure for this long is forgotten.
FORGET_AFTER = 60 * 60.0

#: Pruning bound, so random usernames cannot grow the table without limit.
_MAX_KEYS = 10_000


@dataclass
class _Failures:
    count: int
    last: float
    blocked_until: float


_failures: dict[str, _Failures] = {}

#: Replaced in tests.
_clock: Callable[[], float] = time.monotonic


def _keys(address: str, username: str) -> list[tuple[str, int]]:
    return [
        (f"user:{username.casefold()}", FREE_FAILURES_USER),
        (f"addr:{address}", FREE_FAILURES_ADDRESS),
    ]


def retry_after(address: str, username: str) -> float:
    """Seconds until this address may try this username again (0 = now)."""
    now = _clock()
    wait = 0.0
    for key, _free in _keys(address, username):
        entry = _failures.get(key)
        if entry is not None:
            wait = max(wait, entry.blocked_until - now)
    return max(wait, 0.0)


def record_failure(address: str, username: str) -> None:
    """Count a failed login for both the address and the username."""
    now = _clock()
    if len(_failures) >= _MAX_KEYS:
        _prune(now)
    for key, free in _keys(address, username):
        entry = _failures.get(key)
        if entry is None or now - entry.last > FORGET_AFTER:
            entry = _Failures(count=0, last=now, blocked_until=0.0)
            _failures[key] = entry
        entry.count += 1
        entry.last = now
        if entry.count > free:
            delay = min(BASE_DELAY * 2 ** (entry.count - free - 1), MAX_DELAY)
            entry.blocked_until = now + delay


def record_success(username: str) -> None:
    """Clear the username's failures.

    The address keeps its count: one valid account must not reset the limit
    for guessing the others.
    """
    _failures.pop(f"user:{username.casefold()}", None)


def _prune(now: float) -> None:
    for key in [k for k, e in _failures.items() if now - e.last > FORGET_AFTER]:
        del _failures[key]
    if len(_failures) >= _MAX_KEYS:
        # Still full of recent keys: drop the oldest half.
        for key in sorted(_failures, key=lambda k: _failures[k].last)[: _MAX_KEYS // 2]:
            del _failures[key]


def clear_login_throttle() -> None:
    """Forget every failure (tests)."""
    _failures.clear()
