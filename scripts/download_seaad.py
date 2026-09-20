#!/usr/bin/env python
"""Download and stage the public SEA-AD MTG MERFISH dataset.

Open access (AWS Open Data Registry, no DUA). Downloads the combined
all-donors h5ad (~500 MB) with size + SHA-256 integrity checks and writes
``MANIFEST.json`` into the cache directory, which unblocks
``seaad_niches.io.load_merfish_cells``.

Usage
-----
    python scripts/download_seaad.py --cache-dir data/seaad
    python scripts/download_seaad.py --cache-dir data/seaad --verify-only
    python scripts/download_seaad.py --list
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from seaad_niches.download import (  # noqa: E402
    DEFAULT_H5AD_KEY,
    DEFAULT_H5AD_SHA256,
    DEFAULT_H5AD_SIZE,
    download_merfish_h5ad,
    list_merfish_files,
    verify_manifest,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", default="data/seaad", help="staging directory")
    ap.add_argument("--key", default=DEFAULT_H5AD_KEY, help="S3 key to download")
    ap.add_argument("--sha256", default=DEFAULT_H5AD_SHA256,
                    help="expected SHA-256 (default: recorded value, if pinned)")
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    ap.add_argument("--verify-only", action="store_true",
                    help="only re-verify the cached file against MANIFEST.json")
    ap.add_argument("--list", action="store_true",
                    help="list available files under the MTG bucket prefix")
    args = ap.parse_args()

    if args.list:
        for f in list_merfish_files():
            print(f"{f['size']:>12}  {f['key']}")
        return 0

    if args.verify_only:
        manifest = verify_manifest(args.cache_dir)
        print(f"OK: {manifest['file']} matches MANIFEST.json "
              f"(sha256={manifest['sha256'][:12]}...)")
        return 0

    dest = download_merfish_h5ad(
        args.cache_dir,
        key=args.key,
        expected_size=DEFAULT_H5AD_SIZE if args.key == DEFAULT_H5AD_KEY else None,
        expected_sha256=args.sha256,
        force=args.force,
    )
    print(f"Staged {dest} ({dest.stat().st_size / 1e6:.0f} MB)")
    print(f"Manifest written to {Path(args.cache_dir) / 'MANIFEST.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
