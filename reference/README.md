# Reference implementations

Independent, side-by-side implementations of the Monotone Adjacent Pooling
Algorithm (MAPA), each idiomatic to its language and free of external
dependencies:

- [`python/`](python/) - pure Python (standard library only)
- [`cpp/`](cpp/) - C++17 (no external dependencies)
- [`r/`](r/) - base R only (`testthat` for tests only)
- [`matlab/`](matlab/) - base MATLAB or GNU Octave (free, open source), no toolboxes required
- [`sas/`](sas/) - SAS macros

Python, C++ and R are validated against the same fixture files in
[`fixtures/`](fixtures/) as part of automated CI. The MATLAB and SAS
implementations use the same fixtures but are validated manually (see their
respective READMEs) since those runtimes are not available in CI.

See [`../docs/mapa-methodology.md`](../docs/mapa-methodology.md) for an
explanation of the algorithm.

## Fixtures

- `fixtures/raw_observations.csv` - raw observation-level data:
  `score,bad`, where `bad` is 1 for a default and 0 otherwise.
- `fixtures/expected_initial_bins.csv` - the expected result of grouping
  `raw_observations.csv` into one bin per unique score:
  `score_min,score_max,n_obs,n_bads` (with `score_min == score_max` for
  every row).
- `fixtures/expected_pooled_bins.csv` - the expected result of running MAPA
  on those initial bins with the default (non-increasing bad rate)
  direction: `score_min,score_max,n_obs,n_bads`.
- `fixtures/expected_min_size_bins.csv` - the expected result of applying
  `enforce_minimum_size` to `expected_pooled_bins.csv` with
  `min_obs=50, min_bads=10`: `score_min,score_max,n_obs,n_bads`.
- `fixtures/expected_calibrated_bins.csv` - the expected result of applying
  Bayesian (credibility) adjustment to `expected_min_size_bins.csv` with
  `k=10` and the default prior (overall bad rate):
  `score_min,score_max,n_obs,n_bads,pd`. Note that two adjacent bands in
  this fixture deliberately cross after adjustment - see the "Bayesian
  adjustment" section in
  [`../docs/mapa-methodology.md`](../docs/mapa-methodology.md).
- `fixtures/expected_repooled_calibrated_bins.csv` - the expected result of
  applying `repool_calibrated_bins` to `expected_calibrated_bins.csv`:
  `score_min,score_max,n_obs,n_bads,pd`. This merges the crossing pair from
  `expected_calibrated_bins.csv` back into a single monotone band - see the
  "Re-pooling after shrinkage" section in
  [`../docs/mapa-methodology.md`](../docs/mapa-methodology.md).
- `fixtures/expected_smoothed_pds.csv` - the expected result of applying
  `interpolate_pd` to `expected_repooled_calibrated_bins.csv` for each
  unique score in `raw_observations.csv`: `score,pd`. See the "Smoothing:
  log-odds interpolation" section in
  [`../docs/mapa-methodology.md`](../docs/mapa-methodology.md).
- `fixtures/expected_pooled_bins_confidence.csv` - the expected result of
  running MAPA on `expected_initial_bins.csv` with `min_confidence=0.95`:
  `score_min,score_max,n_obs,n_bads`. See the "Confidence-based pooling"
  section in [`../docs/mapa-methodology.md`](../docs/mapa-methodology.md).

### Weighted fixtures

Parallel set of fixtures for value-weighted observations (see the "Number
vs. value weighting" section in
[`../docs/mapa-methodology.md`](../docs/mapa-methodology.md)):

- `fixtures/raw_observations_weighted.csv` - weighted observation-level
  data: `score,bad,weight`.
- `fixtures/expected_initial_bins_weighted.csv` - expected weighted initial
  bins: `score_min,score_max,n_obs,n_bads,count,count_bads`.
- `fixtures/expected_pooled_bins_weighted.csv` - expected weighted pooled
  bins.
- `fixtures/expected_min_size_bins_weighted.csv` - expected weighted
  minimum-size-enforced bins.
- `fixtures/expected_calibrated_bins_weighted.csv` - expected weighted
  calibrated bins (after Bayesian adjustment).
- `fixtures/expected_repooled_calibrated_bins_weighted.csv` - expected
  weighted re-pooled calibrated bins.
- `fixtures/expected_smoothed_pds_weighted.csv` - expected weighted
  smoothed PDs.
