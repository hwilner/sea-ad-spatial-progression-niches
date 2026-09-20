"""Download and cache the public SEA-AD MTG MERFISH dataset.

The SEA-AD spatial transcriptomics (MERFISH, middle temporal gyrus) data
are *open access* -- no DUA is required. They are distributed by the Allen
Institute on the AWS Open Data Registry bucket ``sea-ad-spatial-transcriptomics``
(registry entry: https://registry.opendata.aws/allen-sea-ad-atlas), as
documented in Gabitto et al., 2024, Nature Neuroscience
(doi:10.1038/s41593-024-01774-5).

This module implements:

- :func:`list_merfish_files` -- list available files in the MTG bucket prefix.
- :func:`download_merfish_h5ad` -- streamed download with SHA-256 integrity
  verification and a ``MANIFEST.json`` cache manifest (the manifest gates
  :mod:`seaad_niches.io` real-data loaders).
- :func:`verify_manifest` -- re-verify cached files against the manifest.

Only the processed, open h5ad is downloaded. Raw/controlled-access material
(AD Knowledge Portal) is out of scope.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

BUCKET_BASE = "https://sea-ad-spatial-transcriptomics.s3.amazonaws.com"
MTG_PREFIX = "middle-temporal-gyrus/"

#: Combined all-donors MERFISH h5ad (processed, open access).
DEFAULT_H5AD_KEY = (
    "middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad"
)
#: Byte size reported by the bucket listing (integrity pre-check).
DEFAULT_H5AD_SIZE = 497_737_587
#: SHA-256 of the default file, verified against a reference download.
#: Recorded in docs/DATA_SOURCES.md; ``None`` disables the hash check.
DEFAULT_H5AD_SHA256: str | None = (
    "3e3dac22446a8ce66afd07c209cd74df260b2f3d7c4085ffe4732390d2054a24"
)

DATASET_CITATION = (
    "Gabitto, Travaglini et al. (2024) Integrated multimodal cell atlas of "
    "Alzheimer's disease. Nature Neuroscience 27, 2366-2383. "
    "doi:10.1038/s41593-024-01774-5"
)

_S3_NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"


def list_merfish_files(prefix: str = MTG_PREFIX, max_keys: int = 1000) -> list[dict]:
    """List files under a bucket prefix (public, unsigned S3 listing).

    Returns a list of dicts with ``key``, ``size`` and ``etag`` (the S3
    ETag; MD5 for single-part uploads).
    """
    url = f"{BUCKET_BASE}/?list-type=2&prefix={prefix}&max-keys={max_keys}"
    with urllib.request.urlopen(url, timeout=120) as resp:
        root = ET.fromstring(resp.read())
    out = []
    for c in root.iter(f"{_S3_NS}Contents"):
        out.append(
            {
                "key": c.findtext(f"{_S3_NS}Key"),
                "size": int(c.findtext(f"{_S3_NS}Size")),
                "etag": c.findtext(f"{_S3_NS}ETag", "").strip('"'),
            }
        )
    return out


def _sha256(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download_merfish_h5ad(
    cache_dir: str | Path,
    key: str = DEFAULT_H5AD_KEY,
    expected_size: int | None = DEFAULT_H5AD_SIZE,
    expected_sha256: str | None = None,
    force: bool = False,
) -> Path:
    """Download the SEA-AD MTG MERFISH h5ad into ``cache_dir`` with checks.

    Streams the file, verifies byte size and (when given) SHA-256, then
    writes ``cache_dir/MANIFEST.json`` recording the dataset key, size,
    checksum, source URL and citation. The manifest is what unblocks the
    real-data branch of :func:`seaad_niches.io.load_merfish_cells`.

    Returns the path to the downloaded h5ad.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / Path(key).name
    url = f"{BUCKET_BASE}/{key}"

    if dest.exists() and not force:
        if expected_size is not None and dest.stat().st_size != expected_size:
            dest.unlink()  # incomplete earlier download; re-fetch
        else:
            _write_manifest(cache_dir, key, dest, url)
            return dest

    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "seaad_niches"})
    with urllib.request.urlopen(req, timeout=300) as resp, open(tmp, "wb") as f:
        while True:
            block = resp.read(1 << 22)
            if not block:
                break
            f.write(block)
    tmp.rename(dest)

    if expected_size is not None and dest.stat().st_size != expected_size:
        raise IOError(
            f"Size mismatch for {key}: got {dest.stat().st_size}, "
            f"expected {expected_size}"
        )
    digest = _sha256(dest)
    if expected_sha256 is not None and digest != expected_sha256:
        raise IOError(
            f"SHA-256 mismatch for {key}: got {digest}, "
            f"expected {expected_sha256}"
        )
    _write_manifest(cache_dir, key, dest, url, sha256=digest)
    return dest


def _write_manifest(
    cache_dir: Path, key: str, dest: Path, url: str, sha256: str | None = None
) -> Path:
    manifest = {
        "dataset": "SEA-AD MTG MERFISH",
        "key": key,
        "source_url": url,
        "registry": "https://registry.opendata.aws/allen-sea-ad-atlas",
        "citation": DATASET_CITATION,
        "file": dest.name,
        "size_bytes": dest.stat().st_size,
        "sha256": sha256 or _sha256(dest),
        "downloaded_unix": int(time.time()),
        "access": "open (AWS Open Data Registry; no DUA required)",
    }
    path = cache_dir / "MANIFEST.json"
    path.write_text(json.dumps(manifest, indent=2))
    return path


def verify_manifest(cache_dir: str | Path) -> dict:
    """Re-hash cached files and compare against MANIFEST.json.

    Returns the manifest dict; raises IOError on any mismatch.
    """
    cache_dir = Path(cache_dir)
    manifest = json.loads((cache_dir / "MANIFEST.json").read_text())
    dest = cache_dir / manifest["file"]
    if not dest.exists():
        raise IOError(f"Cached file missing: {dest}")
    if dest.stat().st_size != manifest["size_bytes"]:
        raise IOError("Cached file size does not match MANIFEST.json")
    if manifest.get("sha256") and _sha256(dest) != manifest["sha256"]:
        raise IOError("Cached file SHA-256 does not match MANIFEST.json")
    return manifest
