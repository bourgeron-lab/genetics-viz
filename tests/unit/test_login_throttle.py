"""Login backoff (utils/login_throttle.py)."""

import pytest

from genetics_viz.utils import login_throttle as throttle


@pytest.fixture
def clock(monkeypatch):
    now = {"t": 1000.0}
    monkeypatch.setattr(throttle, "_clock", lambda: now["t"])
    return now


def _fail(n, address="10.0.0.1", username="alice"):
    for _ in range(n):
        throttle.record_failure(address, username)


def test_free_failures_then_doubling(clock):
    _fail(throttle.FREE_FAILURES_USER)
    assert throttle.retry_after("10.0.0.1", "alice") == 0
    _fail(1)
    assert throttle.retry_after("10.0.0.1", "alice") == throttle.BASE_DELAY
    _fail(1)
    assert throttle.retry_after("10.0.0.1", "alice") == 2 * throttle.BASE_DELAY
    clock["t"] += 2 * throttle.BASE_DELAY
    assert throttle.retry_after("10.0.0.1", "alice") == 0


def test_delay_is_capped(clock):
    _fail(throttle.FREE_FAILURES_USER + 40)
    assert throttle.retry_after("10.0.0.1", "alice") == throttle.MAX_DELAY


def test_username_limit_holds_across_addresses(clock):
    for i in range(throttle.FREE_FAILURES_USER + 1):
        throttle.record_failure(f"10.0.0.{i}", "alice")
    assert throttle.retry_after("10.9.9.9", "ALICE") > 0
    assert throttle.retry_after("10.9.9.9", "bob") == 0


def test_address_limit_holds_across_usernames(clock):
    for i in range(throttle.FREE_FAILURES_ADDRESS + 1):
        throttle.record_failure("10.0.0.1", f"user{i}")
    assert throttle.retry_after("10.0.0.1", "someone-else") > 0
    assert throttle.retry_after("10.0.0.2", "someone-else") == 0


def test_success_clears_the_user_but_not_the_address(clock):
    for i in range(throttle.FREE_FAILURES_ADDRESS + 1):
        throttle.record_failure("10.0.0.1", "alice" if i < 3 else f"u{i}")
    throttle.record_success("alice")
    assert throttle.retry_after("10.0.0.2", "alice") == 0
    assert throttle.retry_after("10.0.0.1", "alice") > 0


def test_old_failures_are_forgotten(clock):
    _fail(throttle.FREE_FAILURES_USER)
    clock["t"] += throttle.FORGET_AFTER + 1
    _fail(1)
    assert throttle.retry_after("10.0.0.1", "alice") == 0


def test_table_is_bounded(clock, monkeypatch):
    monkeypatch.setattr(throttle, "_MAX_KEYS", 10)
    for i in range(50):
        throttle.record_failure("10.0.0.1", f"user{i}")
    assert len(throttle._failures) <= 10
