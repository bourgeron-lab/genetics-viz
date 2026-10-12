"""Multi-data-directory registry (utils/data.py) and auth helpers (utils/auth.py)."""

from types import SimpleNamespace

import pytest

from genetics_viz.config_model import DataDirectoryConfig
from genetics_viz.utils import auth, data


def test_registry_assigns_prefixes_and_default(data_dir, tmp_path):
    other = tmp_path / "other"
    (other / "cohorts").mkdir(parents=True)
    data.init_all_data_stores(
        [
            DataDirectoryConfig(str(other), "second"),
            DataDirectoryConfig(str(data_dir), "main", default=True),
        ]
    )
    assert data.get_static_prefix_map() == {
        str(other): "/data-0",
        str(data_dir): "/data-1",
    }
    assert data.get_default_data_dir_path() == str(data_dir)
    assert "DEMO" in data.get_data_store().cohorts  # no request context: the default
    assert data.get_static_prefix() == "/data-1"
    assert {o["label"] for o in data.get_data_dir_options()} == {
        "other — second",
        "data — main",
    }


def test_empty_registry():
    with pytest.raises(RuntimeError):
        data.get_data_store()
    assert data.get_data_store_or_none() is None
    assert data.get_static_prefix() == "/data"


def test_missing_directory_still_registers(tmp_path):
    data.init_all_data_stores([DataDirectoryConfig(str(tmp_path / "gone"))])
    assert data.get_data_store().cohorts == {}


def test_hot_added_prefixes_never_collide(tmp_path, monkeypatch):
    mounted = []
    import nicegui

    monkeypatch.setattr(
        nicegui.app, "add_static_files", lambda prefix, path: mounted.append(prefix)
    )
    dirs = []
    for name in "abc":
        (tmp_path / name / "cohorts").mkdir(parents=True)
        dirs.append(str(tmp_path / name))
    data.add_data_store(dirs[0])
    data.add_data_store(dirs[1])
    data.remove_data_store(dirs[0])
    data.add_data_store(dirs[2])
    # A static route cannot be unmounted, so a prefix must never be reused.
    assert len(set(mounted)) == len(mounted) == 3
    assert data.add_data_store(dirs[1]) == data.get_static_prefix_map()[dirs[1]]


@pytest.fixture
def session(monkeypatch):
    storage = {}
    monkeypatch.setattr(
        auth, "app", SimpleNamespace(storage=SimpleNamespace(user=storage))
    )
    return storage


@pytest.mark.parametrize(
    "role, write, admin",
    [("reader", False, False), ("curator", True, False), ("administrator", True, True)],
)
def test_roles(session, role, write, admin):
    session.update(authenticated=True, role=role, username="u")
    assert auth.can_write() is write and auth.is_admin() is admin
    assert auth.check_auth() is None and auth.get_current_user() == "u"


def test_unauthenticated_redirects_and_defaults_to_reader(session):
    redirect = auth.check_auth()
    assert redirect is not None and redirect.headers["location"] == "/login"
    assert auth.get_current_role() == "reader" and not auth.can_write()
