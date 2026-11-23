# Results

**Generated 2026-11-23 by executing this repository's own test suite against the
current working tree.**

## Verification run

| Metric | Value |
|---|---|
| Tests passed | **32** |
| Tests failed | 0 |
| Errors | 0 |
| Skipped | 0 |
| Wall clock | 1.73 s |
| Environment | Python 3.12, CPU-only, no GPU, optional `anndata` installed |

Reproduce with:

```bash
PYTHONPATH=$PWD:$PWD/src pytest -q
```

## What this establishes

- The code imports cleanly and **32 tests pass with 0 failures**, 0 skipped.
- The suite now includes synthetic recovery of the pre-registered donor-level
  GAM, an end-to-end smoke test of the one-choice sensitivity harness, and
  fixture tests for the Paper 1 dataset-freeze criteria, QC, and manifest
  machinery.
- That is real computational evidence; it is not a biological finding.

## What this does NOT establish

- **No scientific claim is verified by a passing unit test.** These tests check
  that functions behave as specified -- not that a biological or computational
  hypothesis holds.
- **No benchmark, metric, or comparative result was reproduced.** The unit tests
  are not experiments.
- **Coverage was not measured.** A passing suite can still leave large areas
  untested.
- **Statistical significance was not assessed.** No experiment was run, so there
  is nothing to test.

## Negative and unverified results

Recorded explicitly, so they are not mistaken for successes:

- **No niche-progression association was fit here.** The separately committed
  dataset freeze reports measured metadata/QC counts only; it is not a
  biological association result.
- **The test suite passing is not evidence the research is correct.** It is
  evidence the code is internally consistent.

## Reproducing

```bash
PYTHONPATH=$PWD:$PWD/src pytest -q
```
