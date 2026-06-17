# MAPA - C++ reference implementation

C++17, no external dependencies. See [`mapa.hpp`](mapa.hpp) /
[`mapa.cpp`](mapa.cpp) for the implementation and
[`../../docs/mapa-methodology.md`](../../docs/mapa-methodology.md) for the
algorithm description.

## Usage

### Step by step

```cpp
#include "mapa.hpp"

// Raw (score, bad) observations - bad = 1 for a default, 0 otherwise.
// See "Value-weighted observations" below for the (score, bad, weight) variant.
std::vector<std::pair<double, int>> observations = {
    {400, 1},
    {400, 0},
    {405, 0},
    // ...
};

std::vector<mapa::Bin> pooled = mapa::calibrate(observations);

// Further pool any bands that are too small to support a stable estimate.
std::vector<mapa::Bin> sized = mapa::enforce_minimum_size(pooled, 50, 10);

// Shrink each band's bad rate toward the overall portfolio rate.
std::vector<mapa::CalibratedBin> calibrated = mapa::apply_bayesian_adjustment(sized, 10.0);

// Re-pool to restore monotonicity of pd after shrinkage.
calibrated = mapa::repool_calibrated_bins(calibrated);

// Smooth the step function into a continuous PD curve.
double pd = mapa::interpolate_pd(calibrated, 412.0);
```

### Whole pipeline at once

```cpp
mapa::CalibrationResult result = mapa::run_pipeline(observations, 10.0, 50, 10);

// The band table - the typical deliverable for reporting and governance.
for (const auto& b : result.bands) {
    // b.score_min, b.score_max, b.n_obs, b.n_bads, b.pd
}

// A smoothed, continuous PD for an individual score.
double pd = result.pd_for_score(412.0);
```

`calibrate` groups the raw observations into one bin per unique score
(`bins_from_observations`) and then pools adjacent bins (`mapa`) until the
bad rate is monotone. If you already have pre-aggregated score bands, call
`mapa` directly on a `std::vector<Bin>` with `score_min`/`score_max` set to
each band's range.

Bins carry two sets of counts: `n_obs`/`n_bads` (`double` - weighted sums)
and `count`/`count_bads` (`long` - raw observation counts). For unweighted
observations they are identical.

`enforce_minimum_size` is an optional step that further pools any bands
falling short of the given `min_obs`/`min_bads` thresholds, merging each
into whichever neighbour has the closer bad rate, then re-pooling to
restore monotonicity. With the default thresholds (0, 0) it's a no-op.
The `use_counts` parameter (default `true`) controls whether the
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
`pd_for_score()` (smoothed PDs via `interpolate_pd`). The `use_counts`
parameter is passed through to `enforce_minimum_size`. Use whichever
representation suits the consumer.

### Confidence-based pooling

`mapa` (and therefore `calibrate`, `enforce_minimum_size`, and
`run_pipeline`) accepts an optional `min_confidence` parameter (e.g. `0.95`
for 95%). When given, adjacent bins whose bad rates do not differ at that
confidence level (a two-proportion z-test) are merged as well, even if they
don't violate monotonicity:

```cpp
std::vector<mapa::Bin> pooled = mapa::mapa(initial_bins, /*increasing=*/false, /*min_confidence=*/0.95);
```

This produces fewer, larger bins whose bad rates are more reliably
distinguishable from their neighbours. With the default (`std::nullopt`),
only monotonicity violations are merged, as before.

### Value-weighted observations

Observations can carry an optional weight (typically exposure at default).
Pass a `std::vector<std::tuple<double, int, double>>` instead of the
default `std::vector<std::pair<double, int>>`:

```cpp
// Value-weighted observations: (score, bad, weight)
std::vector<std::tuple<double, int, double>> weighted_obs = {
    {400, 1, 50000.0},
    {400, 0, 120000.0},
    {405, 0, 85000.0},
    // ...
};
mapa::CalibrationResult result = mapa::run_pipeline(weighted_obs, 10.0, 50, 10);
```

When weights are provided, `Bin.n_obs`/`n_bads` are `double` (weighted
sums) and `Bin.count`/`count_bads` are `long` (raw observation counts).
For unweighted observations the two sets are identical. The z-test in
confidence-based pooling always uses `count`/`count_bads` for sample sizes.

## Building and running the tests

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```
