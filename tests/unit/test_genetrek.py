"""GeneTrek snapshot: indexing, identifier resolution, refresh and cache."""

import io
import json
import urllib.request
from pathlib import Path

import datagen
import pytest
from genetics_viz import config_model
from genetics_viz.utils import genetrek


def test_resolution_layers(snapshot):
    assert snapshot.resolve("SHANK3") == "SHANK3"
    assert snapshot.resolve("shank3") == "SHANK3"
    assert snapshot.resolve("ENSG00000251322") == "SHANK3"
    assert snapshot.resolve("ENSG00000251322.12") == "SHANK3"  # versioned ENSG
    assert snapshot.resolve("KIAA1109") == "BLTP1"  # previous symbol
    assert snapshot.resolve("PROSAP2") == "SHANK3"  # unambiguous alias
    assert snapshot.resolve("SHARED") is None  # alias of two genes
    assert snapshot.resolve("NOPE") is None
    assert snapshot.resolve("") is None and snapshot.resolve(None) is None


def test_catalog_membership_and_info(snapshot):
    names = [gs.name for gs in snapshot.catalog]
    assert names == ["hc_ndd_v7", "hc_ndd_v6", "sfari_1", "constraint", "syngo"]
    hc = snapshot.geneset("hc_ndd")  # by alias
    assert hc is snapshot.geneset("hc_ndd_v7") and hc.key == "hc_ndd"
    assert snapshot.geneset("hc_ndd_v6").key == "hc_ndd_v6"
    assert snapshot.members["SHANK3"] >= {"hc_ndd_v7", "constraint", "syngo"}
    assert "PLAIN" not in snapshot.members
    info = snapshot.info["SHANK3"]
    assert dict(info.values)["EAGLE"] == "Definitive"
    assert dict(info.values)["LOEUF"] == "0.36"
    assert [f.name for f in snapshot.inheritance_fields] == [
        "inheritance",
        "inheritance_panelapp_t",
    ]


def test_matching_identifiers_expand_old_symbols(snapshot):
    ids = snapshot.matching_identifiers(snapshot.set_members["hc_ndd_v7"])
    assert {"SHANK3", "BLTP1", "KIAA1109", "PROSAP2"} <= ids
    assert "SHARED" not in ids


def test_links(app_config):
    assert genetrek.gene_url("SHANK3") == "http://genetrek.test/?gene=SHANK3"
    url = genetrek.search_url(["SHANK3", "chr22:1-100"])
    assert url == "http://genetrek.test/#search=SHANK3,chr22:1-100&genome=GRCh38"


def test_url_defaults_and_trailing_slash(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text("storage_secret: s\ngenetrek_url: https://example.org/gt\n")
    config_model.load_config(path)
    assert config_model.get_genetrek_url() == "https://example.org/gt/"
    path.write_text("storage_secret: s\n")
    config_model.load_config(path)
    assert config_model.get_genetrek_url() == config_model.DEFAULT_GENETREK_URL
    assert config_model.get_genetrek_cache_dir() == tmp_path / ".genetrek"


def test_status_before_and_without_data():
    snap = genetrek.get_snapshot()
    assert not snap.ok and snap.loading
    assert "Loading" in snap.status_text()


def test_cached_snapshot_loads(app_config):
    snap = genetrek.load_cached_snapshot()
    assert snap is not None and snap.ok and snap.release == datagen.GENETREK_RELEASE
    assert genetrek.get_snapshot() is snap


def test_refresh_failure_keeps_the_snapshot(app_config):
    genetrek.load_cached_snapshot()
    snap = genetrek.refresh_snapshot()  # urlopen is refused in tests
    assert snap.ok and snap.release == datagen.GENETREK_RELEASE
    assert snap.error and "Last check failed" in snap.status_text()


class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _fake_api(monkeypatch, release: str, calls: list) -> None:
    buf = io.BytesIO()
    datagen.genetrek_frame().write_parquet(buf)
    bodies = {
        "api/health": json.dumps({"status": "ok", "release": release}).encode(),
        "api/fields": json.dumps({"fields": datagen.GENETREK_FIELDS}).encode(),
        "api/download": buf.getvalue(),
    }

    def fake_urlopen(request, timeout=None):
        url = request.full_url
        calls.append(url)
        for key, body in bodies.items():
            if key in url:
                return _Response(body)
        raise AssertionError(url)

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)


def test_refresh_downloads_a_new_release_and_caches_it(app_config, monkeypatch):
    genetrek.load_cached_snapshot()
    calls: list = []
    _fake_api(monkeypatch, "new-release", calls)
    snap = genetrek.refresh_snapshot()
    assert snap.release == "new-release" and snap.ok and snap.error is None
    assert len(calls) == 3
    cache = config_model.get_genetrek_cache_dir()
    assert sorted(p.name for p in cache.iterdir()) == [
        "fields.json",
        "new-release.parquet",
    ]

    calls.clear()
    assert genetrek.refresh_snapshot().release == "new-release"
    assert len(calls) == 1, "an unchanged release costs one health call"


def test_unreadable_cache_is_ignored(app_config):
    cache: Path = config_model.get_genetrek_cache_dir()
    (cache / "fields.json").write_text("{not json")
    assert genetrek.load_cached_snapshot() is None


@pytest.mark.parametrize("ident", ["SHANK3", "KIAA1109"])
def test_resolved_symbols_are_current(snapshot, ident):
    assert snapshot.resolve(ident) in snapshot.by_symbol.values()
