# MAPA — R Reference Implementation

A pure base-R implementation of the Monotone Adjacent Pooling Algorithm (MAPA)
for PD calibration. No external packages are required to run the core
algorithm; `testthat` is used only for the test suite.

## Files

| File | Purpose |
|------|---------|
| `mapa.R` | Core implementation — source this in your own scripts |
| `test_mapa.R` | `testthat` test suite mirroring the Python tests |

Fixture data shared with all reference implementations lives in
`../fixtures/`.

## Running the tests

```bash
cd reference/r
Rscript test_mapa.R
```

`testthat` must be installed (`install.packages("testthat")`).

## Usage

### Step-by-step

```r
source("mapa.R")

# Load raw observations — a data.frame with columns `score` and `bad`
# (and optionally `weight` — see "Value-weighted observations" below)
obs <- read.csv("my_observations.csv")

# 1. Group into one bin per unique score
initial <- bins_from_observations(obs)

# 2. Pool until bad rates are monotone (non-increasing by default)
pooled <- mapa(initial)

# 3. Enforce minimum bin sizes
sized <- enforce_minimum_size(pooled, min_obs = 50, min_bads = 10)

# 4. Apply Bayesian shrinkage toward the overall bad rate
calibrated <- apply_bayesian_adjustment(sized, k = 10)

# 5. Restore monotonicity of the adjusted PDs
repooled <- repool_calibrated_bins(calibrated)

# 6. Smooth into a continuous PD curve
pd_at_550 <- interpolate_pd(repooled, score = 550)
```

### Whole-pipeline shortcut

```r
source("mapa.R")

obs      <- read.csv("my_observations.csv")
pipeline <- run_pipeline(obs, k = 10, min_obs = 50, min_bads = 10)

# Step-function band table (data.frame)
print(pipeline$bands)

# Smoothed PD for an individual score
pipeline$pd_for_score(550)
```

### Value-weighted observations

If the input `data.frame` contains a `weight` column, observations are
value-weighted. Without it, all weights default to 1 (number-weighted)
and `n_obs`/`n_bads` equal `count`/`count_bads`.

```r
# Value-weighted: data.frame with score, bad, weight columns
obs <- read.csv("my_weighted_observations.csv")  # has score, bad, weight
pipeline <- run_pipeline(obs, k = 10, min_obs = 50, min_bads = 10, use_counts = TRUE)
```

The z-test in confidence-based pooling always uses `count`/`count_bads`
for sample sizes.

## Optional parameters

| Parameter | Function(s) | Default | Meaning |
|-----------|-------------|---------|---------|
| `increasing` | `mapa`, `calibrate`, `enforce_minimum_size`, `repool_calibrated_bins`, `run_pipeline` | `FALSE` | `TRUE` makes bad rate non-decreasing (higher score = higher risk) |
| `min_confidence` | `mapa`, `calibrate`, `enforce_minimum_size`, `run_pipeline` | `NULL` | Confidence level (e.g. `0.95`) for the two-proportion z-test; adjacent bins whose rates are not significantly different are merged |
| `prior` | `apply_bayesian_adjustment`, `run_pipeline` | `NULL` (overall bad rate) | PD to shrink toward in Bayesian adjustment |
| `use_counts` | `enforce_minimum_size`, `run_pipeline` | `TRUE` | When `TRUE`, `min_obs`/`min_bads` thresholds check raw observation counts (`count`/`count_bads`); when `FALSE`, they check weighted sums (`n_obs`/`n_bads`) |

## Bin data structure

Bins are plain `data.frame`s with these columns:

| Column | Type | Notes |
|--------|------|-------|
| `score_min` | numeric | Lower bound of the score band |
| `score_max` | numeric | Upper bound of the score band |
| `n_obs` | numeric | Weighted sum of observations in the band |
| `n_bads` | numeric | Weighted sum of defaults (bad == 1) |
| `count` | integer | Raw number of observations in the band |
| `count_bads` | integer | Raw number of defaults (bad == 1) |
| `pd` | numeric | Bayesian-adjusted PD (calibrated bins only) |
