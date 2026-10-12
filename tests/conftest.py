"""Shared fixtures: synthetic data, an application config, and a clean global state.

Most of the app keeps module-level state (the loaded config, the data store
registry, file caches keyed on stat, a gene-scoring singleton). Without the
autouse reset below, tests would pass or fail depending on their order.
"""

import hashlib
from pathlib import Path
from typing import Iterator

import pytest
import yaml

import datagen
from genetics_viz import config_model
from genetics_viz.utils import (
    ancestry,
    change_monitor,
    cohort_state,
    data,
    gene_scoring,
    ghfc_ngs_docs,
    pipeline_params,
    roh,
    str_eh,
    sharding,
)

#: Test accounts: username -> (password, role). Test values only.
USERS = {
    "admin": ("admin-pw", "administrator"),
    "curator": ("curator-pw", "curator"),
    "reader": ("reader-pw", "reader"),
}


def _reset_globals() -> None:
    config_model._app_config = None
    data.init_all_data_stores([])
    ancestry.clear_reference_cache()
    roh.clear_reference_cache()
    str_eh.clear_eh_cache()
    cohort_state.clear_state_cache()
    cohort_state.clear_sha_cache()
    pipeline_params.clear_params_cache()
    ghfc_ngs_docs.clear_doc_cache()
    sharding.clear_sharding_cache()
    gene_scoring.GeneScorer._instance = None
    gene_scoring.GeneScorer._initialized = False
    gene_scoring._gene_scorer = None
    change_monitor._change_callbacks.clear()


@pytest.fixture(autouse=True)
def clean_globals() -> Iterator[None]:
    _reset_globals()
    yield
    _reset_globals()


@pytest.fixture(autouse=True)
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail loudly if anything tries to fetch the ghfc-ngs docs for real."""

    def refuse(*args, **kwargs):
        raise AssertionError("network access attempted in a test")

    monkeypatch.setattr(ghfc_ngs_docs.urllib.request, "urlopen", refuse)


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    """A synthetic data directory (see datagen.py)."""
    return datagen.build_data_dir(tmp_path)


@pytest.fixture
def reference_root(tmp_path: Path) -> Path:
    """A synthetic reference bundle root holding ``v1.1.0/bundle``."""
    return datagen.build_reference(tmp_path / "references")


@pytest.fixture
def config_file(tmp_path: Path, data_dir: Path, reference_root: Path) -> Path:
    """An application YAML config over the synthetic data and reference."""
    path = tmp_path / "genetics_viz.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "data_directories": [
                    {"path": str(data_dir), "description": "synthetic", "default": True}
                ],
                "user_list": [
                    {
                        "username": name,
                        "password": hashlib.sha512(pw.encode()).hexdigest(),
                        "role": role,
                    }
                    for name, (pw, role) in USERS.items()
                ],
                "storage_secret": "test-storage-secret",
                "poll_interval": 3600,
                "ancestry_reference_dir": str(reference_root),
            }
        )
    )
    return path


@pytest.fixture
def app_config(config_file: Path) -> config_model.AppConfig:
    """Load the config into ``config_model``'s singleton, as the app does."""
    return config_model.load_config(config_file)


# ---------------------------------------------------------------------------
# Simulated browser (NiceGUI user simulation)
# ---------------------------------------------------------------------------

UI_MAIN = Path(__file__).parent / "ui" / "app_main.py"


@pytest.fixture
async def user(
    config_file: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog
):
    """A simulated user on a freshly booted app over the synthetic data.

    Fails the test on any ERROR log: NiceGUI reports exceptions raised in event
    handlers and timers that way, so a lazy tab that crashes while loading
    fails here even though no exception reaches the test.
    """
    from nicegui.storage import Storage
    from nicegui.testing import user_simulation

    monkeypatch.setenv("GENETICS_VIZ_TEST_CONFIG", str(config_file))
    # NiceGUI's reset deletes Storage.path, which defaults to the real ./.nicegui.
    monkeypatch.setattr(Storage, "path", tmp_path / ".nicegui")
    async with user_simulation(main_file=UI_MAIN) as simulated:
        yield simulated
    errors = [r for r in caplog.get_records("call") if r.levelname == "ERROR"]
    assert not errors, "\n".join(r.getMessage() for r in errors)
