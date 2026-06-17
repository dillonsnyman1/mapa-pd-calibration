# MAPA - Python reference implementation

Pure Python (standard library only). See
[`mapa.py`](mapa.py) for the implementation and
[`../../docs/mapa-methodology.md`](../../docs/mapa-methodology.md) for the
algorithm description.

## Usage

### Step by step

```python
from mapa import (
    apply_bayesian_adjustment,
    calibrate,
    enforce_minimum_size,
    interpolate_pd,
    repool_calibrated_bins,
)

# Raw (score, bad) observations - bad = 1 for a default, 0 otherwise.
# Each tuple is (score, bad) or (score, bad, weight) — see "Value-weighted
# observations" below.
observations = [
    (400, 1),
    (400, 0),
    (405, 0),
    # ...
]

pooled = calibrate(observations)

for b in pooled:
    print(b.score_min, b.score_max, b.n_obs, b.n_bads, b.count, b.count_bads, b.bad_rate)

# Further pool any bands that are too small to support a stable estimate.
sized = enforce_minimum_size(pooled, min_obs=50, min_bads=10)

# Shrink each band's bad rate toward the overall portfolio rate.
calibrated = apply_bayesian_adjustment(sized, k=10)

# Re-pool to restore monotonicity of pd after shrinkage.
calibrated = repool_calibrated_bins(calibrated)

for b in calibrated:
    print(b.score_min, b.score_max, b.n_obs, b.n_bads, b.pd)

# Smooth the step function into a continuous PD curve.
print(interpolate_pd(calibrated, score=412))
```

### Whole pipeline at once

```python
from mapa import run_pipeline

result = run_pipeline(observations, k=10, min_obs=50, min_bads=10)

# The band table - the typical deliverable for reporting and governance.
for b in result.bands:
    print(b.score_min, b.score_max, b.n_obs, b.n_bads, b.pd)

# A smoothed, continuous PD for an individual score.
print(result.pd_for_score(412))
```

`calibrate` groups the raw observations into one bin per unique score
(`bins_from_observations`) and then pools adjacent bins (`mapa`) until the
bad rate is monotone. If you already have pre-aggregated score bands, call
`mapa` directly on a list of `Bin(score_min, score_max, n_obs, n_bads)`.

Bins carry two sets of counts: `n_obs`/`n_bads` (floats - weighted sums)
and `count`/`count_bads` (ints - raw observation counts). For unweighted
observations they are identical.

`enforce_minimum_size` is an optional step that further pools any bands
falling short of the given `min_obs`/`min_bads` thresholds, merging each
into whichever neighbour has the closer bad rate, then re-pooling to
restore monotonicity. With the default thresholds (0, 0) it's a no-op.
The `use_counts` parameter (default `True`) controls whether the
thresholds are checked against the raw observation counts
(`count`/`count_bads`) or against the weighted sums
(`n_obs`/`n_bads`).

`apply_bayesian_adjustment` is an optional step that shrinks each band's
bad rate toward a prior (by default the overall bad rate) using
credibility weighting. Note that this can re-introduce small monotonicity
violations - see
[`../../docs/mapa-methodology.md`](../../docs/mapa-methodology.md) for
details.

`repool_calibrated_bins` is an optional step that re-applies pooling to the
Bayesian-adjusted bands, this time on `pd`, merging any bands whose `pd`
violates monotonicity by taking an `n_obs`-weighted average of their `pd`
values. This restores the monotonicity that `apply_bayesian_adjustment` can
disturb.

`interpolate_pd` is an optional final step that turns the pooled PD step
function into a continuous curve: each pool is reduced to an anchor point
(its midpoint score and the log-odds of its `pd`), and an individual
score's PD is found by linearly interpolating log-odds between the two
bracketing pools' anchors (flat extrapolation beyond the first/last
anchor). Because log-odds is a monotonic transform of `pd`, this preserves
monotonicity.

`run_pipeline` chains all of the above (`calibrate`, `enforce_minimum_size`,
`apply_bayesian_adjustment`, `repool_calibrated_bins`) and returns a
`CalibrationResult` bundling both outputs: `bands` (the band table) and
`pd_for_score` (smoothed PDs via `interpolate_pd`). The `use_counts`
parameter is passed through to `enforce_minimum_size`. Use whichever
representation suits the consumer.

### Confidence-based pooling

`mapa` (and therefore `calibrate`, `enforce_minimum_size`, and
`run_pipeline`) accepts an optional `min_confidence` parameter (e.g. `0.95`
for 95%). When given, adjacent bins whose bad rates do not differ at that
confidence level (a two-proportion z-test) are merged as well, even if they
don't violate monotonicity:

```python
pooled = mapa(initial_bins, min_confidence=0.95)
```

This produces fewer, larger bins whose bad rates are more reliably
distinguishable from their neighbours. With the default (`None`), only
monotonicity violations are merged, as before.

### Value-weighted observations

Observations can carry an optional weight (typically exposure at default).
Pass each observation as a `(score, bad, weight)` tuple instead of
`(score, bad)`. When no weight is given, it defaults to 1
(number-weighted), and `n_obs`/`n_bads` equal `count`/`count_bads`.

```python
# Value-weighted observations: (score, bad, weight)
# Weight is typically exposure at default (EAD).
observations = [
    (400, 1, 50000),
    (400, 0, 120000),
    (405, 0, 85000),
    # ...
]

result = run_pipeline(observations, k=10, min_obs=50, min_bads=10, use_counts=True)
```

When `use_counts=True` (the default), the `min_obs`/`min_bads` thresholds
in `enforce_minimum_size` are checked against the raw observation counts
(`count`/`count_bads`). Set `use_counts=False` to check against the
weighted sums (`n_obs`/`n_bads`) instead. The z-test in confidence-based
pooling always uses the raw counts (`count`/`count_bads`) for sample sizes.

## Running the tests

```bash
pip install pytest
pytest
```
