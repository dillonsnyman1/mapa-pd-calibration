# MAPA - SAS reference implementation

A SAS macro library implementing the same pipeline as the Python and C++
reference implementations:

```
%mapa_bins_from_observations  -> %mapa_pool  (together: %mapa_calibrate)
                                       |
                                       v
                          %mapa_enforce_minimum_size
                                       |
                                       v
                          %mapa_bayesian_adjustment
                                       |
                                       v
                          %mapa_repool_calibrated  -----------> band table (out_bands)
                                       |
                                       v
                          %mapa_interpolate_pd  --------------> smoothed PDs (out_smoothed)
```

`%mapa_run_pipeline` runs the whole chain (from raw observations through
`%mapa_repool_calibrated`, plus `%mapa_interpolate_pd` if `scores`/
`out_smoothed` are given) in one call, producing both outputs above.

See [`mapa.sas`](mapa.sas) for the macro definitions and
[`../../docs/mapa-methodology.md`](../../docs/mapa-methodology.md) for the
algorithm description.

> **Note on testing**: unlike the Python and C++ implementations, this SAS
> version is **not** run as part of automated CI (no SAS runtime is
> available in CI). It has been validated manually against the shared
> fixtures in [`../fixtures/`](../fixtures/) using
> [`test_mapa.sas`](test_mapa.sas).

## Usage

```sas
%include "mapa.sas";

/* raw_obs has columns: score, bad (1 = default, 0 otherwise) */
%mapa_calibrate(in=raw_obs, out=pooled_bins)

%mapa_enforce_minimum_size(in=pooled_bins, out=sized_bins, min_obs=50, min_bads=10)

%mapa_bayesian_adjustment(in=sized_bins, out=calibrated_bins, k=10)

%mapa_repool_calibrated(in=calibrated_bins, out=repooled_bins)

/* scores has one column: score */
%mapa_interpolate_pd(bins=repooled_bins, scores=scores, out=smoothed_pds)
```

Each macro reads/writes datasets with columns `score_min`, `score_max`,
`n_obs`, `n_bads` (plus `pd` for `%mapa_bayesian_adjustment`'s and
`%mapa_repool_calibrated`'s output). If you already have pre-aggregated
score bands, skip `%mapa_bins_from_observations` / `%mapa_calibrate` and
call `%mapa_pool` directly on a dataset with those columns.

`%mapa_interpolate_pd` is the odd one out: it doesn't operate on pools, but
on individual scores. `scores` is a dataset with a single column, `score`,
one row per score to compute a smoothed PD for; `smoothed_pds` gets columns
`score` and `pd`.

### Confidence-based pooling

`%mapa_pool` (and therefore `%mapa_calibrate`, `%mapa_enforce_minimum_size`
and `%mapa_run_pipeline`) accepts an optional `min_confidence=` parameter
(e.g. `0.95` for 95%). When given, adjacent bins whose bad rates do not
differ at that confidence level (a two-proportion z-test) are merged as
well, even if they don't violate monotonicity:

```sas
%mapa_pool(in=initial_bins, out=pooled_bins, min_confidence=0.95)
```

This produces fewer, larger bins whose bad rates are more reliably
distinguishable from their neighbours. With the default (`min_confidence=0`),
only monotonicity violations are merged, as before.

### Running the whole pipeline in one call

`%mapa_run_pipeline` chains all of the above:

```sas
%mapa_run_pipeline(in=raw_obs, out_bands=bands, k=10, min_obs=50, min_bads=10,
                    scores=scores, out_smoothed=smoothed_pds)
```

This produces two independent outputs - use whichever suits the consumer:

- `bands`: the band table (`score_min`, `score_max`, `n_obs`, `n_bads`,
  `pd`) - the typical deliverable for reporting and governance.
- `smoothed_pds`: a smoothed, continuous PD per score (`score`, `pd`).
  Only produced if both `scores` and `out_smoothed` are given.

## Running the driver against the fixtures

[`test_mapa.sas`](test_mapa.sas) runs the full pipeline against
[`../fixtures/raw_observations.csv`](../fixtures/raw_observations.csv) and
prints each step's output for comparison against the corresponding
`expected_*.csv` fixture.

1. Open `test_mapa.sas` in SAS (e.g. SAS Studio / SAS OnDemand for
   Academics, or a local SAS installation).
2. Update the `fixtures_path` macro variable at the top to the absolute
   path of [`../fixtures`](../fixtures) on your system.
3. Run the script and compare the printed tables to the `expected_*.csv`
   files.
