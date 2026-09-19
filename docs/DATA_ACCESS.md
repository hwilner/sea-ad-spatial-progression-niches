# SEA-AD data access

**Status: blocked.** Real-data tasks in this repository require access to
the SEA-AD atlas (snRNA-seq and MERFISH, middle temporal gyrus). The
repository owner must first complete the data-access/data-use steps
tracked in the backlog issue **"Set up SEA-AD open-data access and
environment"** (task-key `setup-data-access-env`, issue #5).

Until access is in place:

- `seaad_niches.io.load_merfish_cells` and `load_snrnaseq_reference`
  raise `RuntimeError` with a pointer to this page.
- All development and CI run on synthetic tissues from
  `seaad_niches.simulate.simulate_tissue`.
- No SEA-AD data or derived real-data results may be committed.

## Planned sources (open, processed data only)

- SEA-AD snRNA-seq and MERFISH via the CELLxGENE collection.
- AWS Open Data Registry: `allen-sea-ad-atlas`.

Exact dataset versions/DOIs will be recorded in `docs/DATA_SOURCES.md`
as part of the data-access task. Controlled-access raw material is out of
scope for this repository.

## Once access is granted

1. Complete issue #5 (document versions, scripted open-data download).
2. Place a `MANIFEST.json` in the local data cache directory.
3. The `io` module is then implemented against the frozen manifest
   (freeze-analytical-dataset task) — the analysis modules
   (`spatial`, `niches`, `progression`) already accept the coordinate /
   label / donor-score tables the loaders will produce.
