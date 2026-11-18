# Results

**Generated 2026-11-18 by executing this repository's own test suite against the
current `main`.**

## Verification run

| Metric | Value |
|---|---|
| Tests passed | **22** |
| Tests failed | 0 |
| Errors | 0 |
| Skipped | 1 |
| Wall clock | 1.5 s |
| Environment | Python 3.12, CPU-only, no GPU |

Reproduce with:

```bash
PYTHONPATH=$PWD:$PWD/src pytest -q
```

## What this establishes

- The code imports cleanly and **22 tests pass with 0 failures**, 1 skipped.
- The suite now includes synthetic recovery of the pre-registered donor-level
  GAM and an end-to-end smoke test of the one-choice sensitivity harness.
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

- **No empirical result was reproduced for this repository.** Any finding
  described in the README or `docs/` was not re-derived here.
- **The test suite passing is not evidence the research is correct.** It is
  evidence the code is internally consistent.

## Reproducing

```bash
PYTHONPATH=$PWD:$PWD/src pytest -q
```
