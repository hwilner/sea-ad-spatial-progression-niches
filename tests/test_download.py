"""Tests for the download/cache machinery and the manifest-gated io loader.

All network access is mocked or replaced with tiny local fixtures; real
SEA-AD downloads are optional and skipped by default.
"""

import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from seaad_niches import download, io


def _make_fake_h5ad(path, n=30):
    ad = pytest.importorskip("anndata")
    obs = pd.DataFrame(
        {
            "donor": ["D1"] * n,
            "section": ["S1"] * (n // 2) + ["S2"] * (n - n // 2),
            "x": np.linspace(0, 100, n),
            "y": np.linspace(0, 50, n),
            "subclass": (["Astro", "Micro", "Exc", "Inh"] * (n // 4 + 1))[:n],
        },
        index=[f"c{i}" for i in range(n)],
    )
    ad.AnnData(obs=obs).write_h5ad(path)


def test_download_verifies_size_and_writes_manifest(tmp_path, monkeypatch):
    payload = b"fake-h5ad-bytes" * 100
    digest = hashlib.sha256(payload).hexdigest()

    class _Resp:
        def __init__(self, data):
            self._buf = [data]

        def read(self, n=-1):
            return self._buf.pop(0) if self._buf else b""

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(
        download.urllib.request, "urlopen", lambda *a, **k: _Resp(payload)
    )
    dest = download.download_merfish_h5ad(
        tmp_path, key="k/file.h5ad", expected_size=len(payload),
        expected_sha256=digest,
    )
    manifest = json.loads((tmp_path / "MANIFEST.json").read_text())
    assert manifest["sha256"] == digest
    assert manifest["size_bytes"] == len(payload)
    # Re-running without force reuses the cache and verifies fine.
    assert download.download_merfish_h5ad(
        tmp_path, key="k/file.h5ad", expected_size=len(payload)
    ) == dest
    assert download.verify_manifest(tmp_path)["file"] == dest.name


def test_download_rejects_corrupt_payload(tmp_path, monkeypatch):
    payload = b"abc"

    class _Resp:
        def read(self, n=-1):
            b, self._done = getattr(self, "_done", payload), b""
            return b

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(
        download.urllib.request, "urlopen", lambda *a, **k: _Resp()
    )
    with pytest.raises(IOError, match="Size mismatch"):
        download.download_merfish_h5ad(
            tmp_path, key="k/file.h5ad", expected_size=999
        )


def test_verify_manifest_detects_tampering(tmp_path):
    f = tmp_path / "f.bin"
    f.write_bytes(b"0123456789")
    (tmp_path / "MANIFEST.json").write_text(
        json.dumps({"file": "f.bin", "size_bytes": 10,
                    "sha256": hashlib.sha256(b"different!").hexdigest()})
    )
    with pytest.raises(IOError, match="SHA-256"):
        download.verify_manifest(tmp_path)


def test_load_merfish_cells_from_staged_cache(tmp_path):
    h5ad = tmp_path / "mini.h5ad"
    _make_fake_h5ad(h5ad)
    (tmp_path / "MANIFEST.json").write_text(json.dumps({"file": "mini.h5ad"}))
    cells = io.load_merfish_cells(tmp_path)
    assert set(cells.columns) == {"cell_id", "donor", "section", "x", "y", "cell_type"}
    assert cells["section"].nunique() == 2
    assert cells["cell_type"].nunique() == 4
    assert np.isfinite(cells[["x", "y"]].to_numpy()).all()


def test_load_merfish_still_gated_without_manifest(tmp_path):
    with pytest.raises(RuntimeError, match="SEA-AD data access"):
        io.load_merfish_cells(cache_dir=tmp_path)
