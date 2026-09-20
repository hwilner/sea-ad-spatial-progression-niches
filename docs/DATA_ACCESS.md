# SEA-AD data access

**Status: unblocked for the open MERFISH dataset.** The SEA-AD MTG MERFISH
spatial transcriptomics data are open access (AWS Open Data Registry,
`sea-ad-spatial-transcriptomics` bucket) — no DUA is required for the
processed data this repository uses. See
[docs/DATA_SOURCES.md](DATA_SOURCES.md) for exact dataset versions, DOIs,
and integrity checksums.

## One-command setup

```bash
make real-pipeline   # download (if needed) + neighbor enrichment + niche detection
```

or step by step:

```bash
python scripts/download_seaad.py --cache-dir data/seaad   # ~500 MB, checksummed
python scripts/analyze_seaad.py --cache-dir data/seaad    # writes reports/
```

The download script verifies byte size and records SHA-256 in
`data/seaad/MANIFEST.json`. The presence of that manifest unblocks
`seaad_niches.io.load_merfish_cells`, which returns a per-cell table with
`cell_id, donor, section, x, y, cell_type`.

## What is gated vs. open

| Data | Access | Used here |
|---|---|---|
| SEA-AD MTG MERFISH (processed h5ad) | Open (AWS Open Data Registry) | Yes — downloaded by `scripts/download_seaad.py` |
| SEA-AD MTG snRNA-seq (processed) | Open (AWS / CELLxGENE) | Not needed yet (MERFISH cells carry mapped subclass/supertype labels) |
| Raw sequencing data (AD Knowledge Portal) | Controlled access (DUA via Sage) | No — out of scope |

## Rules

- Downloaded data live in `data/seaad/` and must **never** be committed
  (`.gitignore` covers `data/`).
- Only small derived aggregate outputs (CSV/JSON under `reports/`) are
  committed.
- Without a staged `MANIFEST.json`, the `io` loaders raise `RuntimeError`
  and the test suite runs on synthetic tissues only.
