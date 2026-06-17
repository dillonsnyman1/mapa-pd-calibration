/******************************************************************************
 Driver script for the MAPA SAS macros (mapa.sas).

 Runs the full pipeline against the shared fixtures in ../fixtures and
 prints each step's output so it can be compared by hand against:

   - ../fixtures/expected_initial_bins.csv
   - ../fixtures/expected_pooled_bins.csv
   - ../fixtures/expected_min_size_bins.csv             (min_obs=50, min_bads=10)
   - ../fixtures/expected_calibrated_bins.csv           (k=10, applied to sized bins)
   - ../fixtures/expected_repooled_calibrated_bins.csv  (k=10)
   - ../fixtures/expected_smoothed_pds.csv
   - ../fixtures/expected_pooled_bins_confidence.csv     (min_confidence=0.95)

 Value-weighted pipeline (uses raw_observations_weighted.csv):
   - ../fixtures/expected_initial_bins_weighted.csv
   - ../fixtures/expected_pooled_bins_weighted.csv
   - ../fixtures/expected_min_size_bins_weighted.csv
   - ../fixtures/expected_repooled_calibrated_bins_weighted.csv
   - ../fixtures/expected_smoothed_pds_weighted.csv

 Update FIXTURES_PATH below to the absolute path of the `reference/fixtures`
 directory on whatever system you're running this on (e.g. SAS OnDemand for
 Academics).
******************************************************************************/

%let fixtures_path = /path/to/mapa-pd-calibration/reference/fixtures;

%include "mapa.sas";

filename raw "&fixtures_path./raw_observations.csv";

proc import datafile=raw out=raw_obs dbms=csv replace;
    getnames=yes;
run;

/* Step 1: bin raw observations by unique score */
%mapa_bins_from_observations(in=raw_obs, out=initial_bins)

/* Step 2: pool adjacent bins until the bad rate is monotone */
%mapa_pool(in=initial_bins, out=pooled_bins)

/* Step 3: further pool bins below the minimum size thresholds */
%mapa_enforce_minimum_size(in=pooled_bins, out=sized_bins, min_obs=50, min_bads=10)

/* Step 4: shrink each bin's bad rate toward the overall bad rate */
%mapa_bayesian_adjustment(in=sized_bins, out=calibrated_bins, k=10)

/* Step 5: re-pool the calibrated bins to restore monotonicity of pd */
%mapa_repool_calibrated(in=calibrated_bins, out=repooled_calibrated_bins)

/* Step 6: smooth the pooled PD step function via log-odds interpolation */
proc sql;
    create table scores as
    select distinct score as score from raw_obs order by score;
quit;

%mapa_interpolate_pd(bins=repooled_calibrated_bins, scores=scores, out=smoothed_pds)

/* Steps 2-6 in one call: %mapa_run_pipeline should reproduce the same
   band table and smoothed PDs as the individual macro calls above. */
%mapa_run_pipeline(in=raw_obs, out_bands=pipeline_bands, k=10, min_obs=50, min_bads=10,
                    scores=scores, out_smoothed=pipeline_smoothed_pds)

/* Step 7: confidence-based pooling - merge adjacent bins whose bad rates
   aren't statistically distinguishable at the 95% confidence level, even
   if they don't violate monotonicity. */
%mapa_pool(in=initial_bins, out=pooled_bins_confidence, min_confidence=0.95)

title "Initial bins (one per unique score) - compare to expected_initial_bins.csv";
proc print data=initial_bins noobs; run;

title "Pooled bins (MAPA) - compare to expected_pooled_bins.csv";
proc print data=pooled_bins noobs; run;

title "Pooled bins after enforce_minimum_size (min_obs=50, min_bads=10) - compare to expected_min_size_bins.csv";
proc print data=sized_bins noobs; run;

title "Calibrated bins (Bayesian adjustment, k=10) - compare to expected_calibrated_bins.csv";
proc print data=calibrated_bins noobs; run;

title "Calibrated bins after repool_calibrated - compare to expected_repooled_calibrated_bins.csv";
proc print data=repooled_calibrated_bins noobs; run;

title "Smoothed PDs (log-odds interpolation) - compare to expected_smoothed_pds.csv";
proc print data=smoothed_pds noobs; run;

title "Pipeline bands (mapa_run_pipeline) - compare to expected_repooled_calibrated_bins.csv";
proc print data=pipeline_bands noobs; run;

title "Pipeline smoothed PDs (mapa_run_pipeline) - compare to expected_smoothed_pds.csv";
proc print data=pipeline_smoothed_pds noobs; run;

title "Pooled bins with confidence-based pooling (min_confidence=0.95) - compare to expected_pooled_bins_confidence.csv";
proc print data=pooled_bins_confidence noobs; run;

/* =========================================================================
   Value-weighted pipeline
   ========================================================================= */

filename raw_w "&fixtures_path./raw_observations_weighted.csv";

proc import datafile=raw_w out=raw_obs_weighted dbms=csv replace;
    getnames=yes;
run;

/* Weighted Step 1: bin raw observations (weight column detected automatically) */
%mapa_bins_from_observations(in=raw_obs_weighted, out=initial_bins_w)

/* Weighted Step 2: pool */
%mapa_pool(in=initial_bins_w, out=pooled_bins_w)

/* Weighted Step 3: enforce minimum size (use_counts=1, thresholds on raw counts) */
%mapa_enforce_minimum_size(in=pooled_bins_w, out=sized_bins_w, min_obs=50, min_bads=10, use_counts=1)

/* Weighted Step 4: Bayesian adjustment */
%mapa_bayesian_adjustment(in=sized_bins_w, out=calibrated_bins_w, k=10)

/* Weighted Step 5: re-pool */
%mapa_repool_calibrated(in=calibrated_bins_w, out=repooled_calibrated_bins_w)

/* Weighted Step 6: smoothed PDs */
proc sql;
    create table scores_w as
    select distinct score as score from raw_obs_weighted order by score;
quit;

%mapa_interpolate_pd(bins=repooled_calibrated_bins_w, scores=scores_w, out=smoothed_pds_w)

/* Weighted run_pipeline in one call */
%mapa_run_pipeline(in=raw_obs_weighted, out_bands=pipeline_bands_w, k=10, min_obs=50, min_bads=10,
                    scores=scores_w, out_smoothed=pipeline_smoothed_pds_w, use_counts=1)

title "Weighted initial bins - compare to expected_initial_bins_weighted.csv";
proc print data=initial_bins_w noobs; run;

title "Weighted pooled bins - compare to expected_pooled_bins_weighted.csv";
proc print data=pooled_bins_w noobs; run;

title "Weighted sized bins (min_obs=50, min_bads=10, use_counts=1) - compare to expected_min_size_bins_weighted.csv";
proc print data=sized_bins_w noobs; run;

title "Weighted repooled calibrated bins - compare to expected_repooled_calibrated_bins_weighted.csv";
proc print data=repooled_calibrated_bins_w noobs; run;

title "Weighted smoothed PDs - compare to expected_smoothed_pds_weighted.csv";
proc print data=smoothed_pds_w noobs; run;

title "Weighted pipeline bands - compare to expected_repooled_calibrated_bins_weighted.csv";
proc print data=pipeline_bands_w noobs; run;

title "Weighted pipeline smoothed PDs - compare to expected_smoothed_pds_weighted.csv";
proc print data=pipeline_smoothed_pds_w noobs; run;

title;
