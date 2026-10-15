"""Per-user gene badge preferences: fallback chain, guarded writes."""

import pytest
import yaml
from genetics_viz import config_model
from genetics_viz.utils import auth, user_preferences as up


@pytest.fixture
def as_user(monkeypatch):
    """Pretend ``name`` is logged in with ``role``."""

    def login(name: str, role: str = "reader") -> None:
        monkeypatch.setattr(up, "get_current_user", lambda: name)
        monkeypatch.setattr(up, "is_admin", lambda: role == "administrator")

    return login


def _prefs(sets, colors=None) -> up.GenePrefs:
    return up.GenePrefs.from_dict({"sets": sets, "tier_colors": colors})


def test_packaged_default_until_anything_is_saved(app_config):
    prefs = up.get_gene_prefs("reader")
    assert prefs.source == "packaged default"
    assert prefs.sets["hc_ndd"] == 1 and prefs.sets["autism_combined_v7"] == "show"
    assert len(prefs.tier_colors) == up.N_TIERS


def test_user_then_admin_default_then_packaged(app_config, as_user):
    as_user("admin", "administrator")
    up.save_default_gene_prefs(_prefs({"constraint": 1}))
    assert up.get_gene_prefs("reader").source == "admin default"
    assert up.get_gene_prefs("reader").sets == {"constraint": 1}

    as_user("reader")
    up.save_gene_prefs("reader", _prefs({"syngo": 2, "hc_ndd": "show"}))
    assert up.get_gene_prefs("reader").sets == {"syngo": 2, "hc_ndd": "show"}
    assert up.get_gene_prefs("curator").source == "admin default"

    up.reset_gene_prefs("reader")
    assert up.get_gene_prefs("reader").source == "admin default"


def test_writes_are_guarded(app_config, as_user):
    as_user("reader")
    with pytest.raises(PermissionError):
        up.save_gene_prefs("curator", _prefs({}))
    with pytest.raises(PermissionError):
        up.save_default_gene_prefs(_prefs({}))
    with pytest.raises(PermissionError):
        up.remove_user_preferences("curator")


def test_saves_merge_and_keep_other_users(app_config, as_user):
    for name in ("reader", "curator"):
        as_user(name)
        up.save_gene_prefs(name, _prefs({"syngo": 3}))
    as_user("admin", "administrator")
    up.remove_user_preferences("reader")
    raw = yaml.safe_load(config_model.get_preferences_path().read_text())
    assert set(raw["users"]) == {"curator"}
    assert config_model.get_preferences_path().name == "genetics_viz.preferences.yaml"


def test_garbage_is_cleaned_on_read(app_config):
    prefs = up.GenePrefs.from_dict(
        {"sets": {"a": 9, "b": "show", "c": "x", "d": "2"}, "tier_colors": ["#000000"]}
    )
    assert prefs.sets == {"b": "show", "d": 2}
    assert prefs.tier_colors[0] == "#000000" and len(prefs.tier_colors) == 4


def test_unreadable_file_falls_back(app_config):
    config_model.get_preferences_path().write_text("users: [unclosed")
    up.clear_preferences_cache()
    assert up.get_gene_prefs("reader").source == "packaged default"


def test_auth_module_is_the_one_patched():
    # The guards call these names through the module, so the fixture works.
    assert up.get_current_user is auth.get_current_user
