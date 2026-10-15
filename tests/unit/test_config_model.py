"""Application config load/save (config_model.py)."""

import yaml

from genetics_viz import config_model
from genetics_viz.config_model import (
    get_ancestry_reference_dir,
    get_default_data_dir,
    hash_password,
    load_config,
    save_config,
    verify_password,
)


def test_password_round_trip():
    hashed = hash_password("s3cret")
    assert verify_password("s3cret", hashed) and not verify_password("other", hashed)


def test_load_sets_the_singleton(app_config, reference_root):
    assert config_model.get_config() is app_config
    assert get_ancestry_reference_dir() == reference_root
    assert app_config.poll_interval == 3600
    assert {u.role for u in app_config.user_list} == {
        "administrator",
        "curator",
        "reader",
    }


def test_missing_storage_secret_is_generated_once(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text("data_directories: []\n")
    first = load_config(path).storage_secret
    assert first and yaml.safe_load(path.read_text())["storage_secret"] == first
    assert load_config(path).storage_secret == first


def test_empty_keys_load_as_defaults(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text(
        "data_directories:\nuser_list:\nstorage_secret: s\nancestry_reference_dir:\npoll_interval:\n"
    )
    config = load_config(path)
    assert config.data_directories == [] and config.user_list == []
    assert config.ancestry_reference_dir == "" and get_ancestry_reference_dir() is None
    assert config.poll_interval == 30


def test_save_omits_defaults(app_config):
    app_config.poll_interval = 30
    app_config.ancestry_reference_dir = ""
    save_config(app_config)
    raw = yaml.safe_load(app_config.config_path.read_text())
    assert "poll_interval" not in raw and "ancestry_reference_dir" not in raw
    assert load_config(app_config.config_path).user_list == app_config.user_list


def test_default_data_dir_falls_back_to_first(app_config):
    assert get_default_data_dir(app_config) == app_config.data_directories[0].path
    app_config.data_directories[0].default = False
    assert get_default_data_dir(app_config) == app_config.data_directories[0].path


def test_genetrek_keys_survive_an_admin_save(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text(
        "storage_secret: s\ngenetrek_url: https://gt.example/\n"
        "genetrek_cache_dir: /tmp/gt\n"
    )
    config = config_model.load_config(path)
    config_model.save_config(config)
    reloaded = config_model.load_config(path)
    assert reloaded.genetrek_url == "https://gt.example/"
    assert reloaded.genetrek_cache_dir == "/tmp/gt"
