# MAPA — MATLAB Reference Implementation

A pure base-MATLAB implementation of the Monotone Adjacent Pooling Algorithm
(MAPA) for PD calibration. No toolboxes are required — the implementation
uses only `erfinv` (available in base MATLAB since R2012a) in place of
`norminv` from the Statistics and Machine Learning Toolbox.

## Files

```
reference/matlab/
├── bins_from_observations.m      Group raw observations into per-score bins
├── mapa_pool.m                   PAVA-style pooling (= Python's mapa())
├── mapa_calibrate.m              bins_from_observations + mapa_pool
├── enforce_minimum_size.m        Pool bins below size thresholds, re-pool
├── apply_bayesian_adjustment.m   Bayesian shrinkage toward a prior PD
├── repool_calibrated_bins.m      Restore monotonicity of pd after shrinkage
├── interpolate_pd.m              Smooth step-function to continuous PD curve
├── run_pipeline.m                Full pipeline (steps above chained together)
├── test_mapa.m                   Script-based test suite
└── private/
    ├── bin_violates.m            Bad-rate monotonicity check
    ├── pd_violates.m             PD monotonicity check
    ├── not_significant.m         Two-proportion z-test helper
    ├── merge_bins.m              Combine two adjacent bins
    └── merge_calibrated.m        Combine two adjacent calibrated bins
```

Fixture data shared with all reference implementations lives in
`../fixtures/`.

## Running the tests

Add `reference/matlab/` to your MATLAB path, then run:

```matlab
run('test_mapa.m')
```

or from inside the `reference/matlab/` directory:

```matlab
test_mapa
```

The script prints each test result and finishes with `All tests passed.`

## Usage

### Step-by-step

```matlab
% Add the reference/matlab folder to the path once
addpath('/path/to/reference/matlab');

% Load raw observations — a table with columns `score` and `bad`
% (and optionally `weight` — see "Value-weighted observations" below)
obs = readtable('my_observations.csv');

% 1. Group into one bin per unique score
initial = bins_from_observations(obs);

% 2. Pool until bad rates are monotone (non-increasing by default)
pooled = mapa_pool(initial);

% 3. Enforce minimum bin sizes
sized = enforce_minimum_size(pooled, 50, 10);

% 4. Apply Bayesian shrinkage toward the overall bad rate
calibrated = apply_bayesian_adjustment(sized, 10);

% 5. Restore monotonicity of the adjusted PDs
repooled = repool_calibrated_bins(calibrated);

% 6. Smooth into a continuous PD curve
pd_at_550 = interpolate_pd(repooled, 550);
```

### Whole-pipeline shortcut

```matlab
addpath('/path/to/reference/matlab');

obs      = readtable('my_observations.csv');
pipeline = run_pipeline(obs, 10, 50, 10);

% Step-function band table
disp(pipeline.bands)

% Smoothed PD for an individual score
pipeline.pd_for_score(550)
```

### Value-weighted observations

If the input table (or N-by-3 matrix) contains a `weight` column,
observations are value-weighted. Without it, all weights default to 1
(number-weighted) and `n_obs`/`n_bads` equal `count`/`count_bads`.

```matlab
% Value-weighted: table or N-by-3 matrix [score, bad, weight]
obs = readtable('my_weighted_observations.csv');  % has score, bad, weight
pipeline = run_pipeline(obs, 10, 50, 10, [], false, [], true);  % use_counts=true
```

The z-test in confidence-based pooling always uses `count`/`count_bads`
for sample sizes.

## Optional parameters

| Parameter | Function(s) | Default | Meaning |
|-----------|-------------|---------|---------|
| `increasing` | `mapa_pool`, `mapa_calibrate`, `enforce_minimum_size`, `repool_calibrated_bins`, `run_pipeline` | `false` | `true` makes bad rate non-decreasing (higher score = higher risk) |
| `min_confidence` | `mapa_pool`, `mapa_calibrate`, `enforce_minimum_size`, `run_pipeline` | `[]` (disabled) | Confidence level (e.g. `0.95`) for the two-proportion z-test; adjacent bins whose rates are not significantly different are merged |
| `prior` | `apply_bayesian_adjustment`, `run_pipeline` | `[]` (overall bad rate) | PD to shrink toward in Bayesian adjustment |
| `use_counts` | `enforce_minimum_size`, `run_pipeline` | `true` | When `true`, `min_obs`/`min_bads` thresholds check raw observation counts (`count`/`count_bads`); when `false`, they check weighted sums (`n_obs`/`n_bads`) |

All optional arguments can be omitted or passed as `[]` to use their defaults.

## Bin data structure

Bins are MATLAB `table` objects with these variables:

| Variable | Type | Notes |
|----------|------|-------|
| `score_min` | double | Lower bound of the score band |
| `score_max` | double | Upper bound of the score band |
| `n_obs` | double | Weighted sum of observations in the band |
| `n_bads` | double | Weighted sum of defaults (bad == 1) |
| `count` | double | Raw number of observations in the band |
| `count_bads` | double | Raw number of defaults (bad == 1) |
| `pd` | double | Bayesian-adjusted PD (calibrated bins only) |

## Naming note

The main pooling function is named `mapa_pool` (not `mapa`) to avoid a
naming conflict with the `reference/matlab/` folder itself, which MATLAB
would otherwise shadow as a package. `mapa_calibrate` similarly corresponds
to Python's `calibrate`.
