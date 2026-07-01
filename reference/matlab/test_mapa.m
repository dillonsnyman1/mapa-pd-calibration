% TEST_MAPA  Script-based test suite for the MATLAB MAPA implementation.
%
% Run with:
%   cd reference/matlab
%   run('test_mapa.m')
% or from any directory that has reference/matlab on the path:
%   run('/path/to/reference/matlab/test_mapa.m')
%
% Prints "All tests passed." on success. Calls error() on failure.
%
% Local helper functions are in the private/ subdirectory.

% Load Octave's datatypes package for table support (no-op in MATLAB).
if exist('OCTAVE_VERSION', 'builtin')
    pkg load datatypes
end

% Locate fixtures relative to this script's directory.
this_dir     = fileparts(mfilename('fullpath'));
fixtures_dir = fullfile(this_dir, '..', 'fixtures');

BAYESIAN_K     = 10;
MIN_OBS        = 50;
MIN_BADS       = 10;
MIN_CONFIDENCE = 0.95;

% ---------------------------------------------------------------------------
% Test 1: bins_from_observations matches expected
% ---------------------------------------------------------------------------
fprintf('Test 1: bins_from_observations matches expected ... ');
obs      = load_observations(fixtures_dir);
result   = bins_from_observations(obs);
expected = load_bins(fixtures_dir, 'expected_initial_bins.csv');
assert(bins_equal(result, expected), 'bins_from_observations does not match expected');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 2: mapa_calibrate matches expected pooled bins
% ---------------------------------------------------------------------------
fprintf('Test 2: mapa_calibrate matches expected pooled bins ... ');
obs      = load_observations(fixtures_dir);
result   = mapa_calibrate(obs);
expected = load_bins(fixtures_dir, 'expected_pooled_bins.csv');
assert(bins_equal(result, expected), 'mapa_calibrate does not match expected_pooled_bins');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 3: result is monotone non-increasing
% ---------------------------------------------------------------------------
fprintf('Test 3: result is monotone non-increasing ... ');
obs    = load_observations(fixtures_dir);
result = mapa_calibrate(obs);
rates  = result.n_bads ./ result.n_obs;
assert(all(diff(rates) <= 0), 'Bad rates are not non-increasing');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 4: pooling preserves totals
% ---------------------------------------------------------------------------
fprintf('Test 4: pooling preserves totals ... ');
obs     = load_observations(fixtures_dir);
initial = bins_from_observations(obs);
result  = mapa_pool(initial);
assert(sum(result.n_obs)  == sum(initial.n_obs),  'n_obs total mismatch');
assert(sum(result.n_bads) == sum(initial.n_bads), 'n_bads total mismatch');
assert(sum(result.n_obs)  == height(obs),         'n_obs total != nrow(obs)');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 5: enforce_minimum_size matches expected
% ---------------------------------------------------------------------------
fprintf('Test 5: enforce_minimum_size matches expected ... ');
obs    = load_observations(fixtures_dir);
pooled = mapa_calibrate(obs);
result = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
expected = load_bins(fixtures_dir, 'expected_min_size_bins.csv');
assert(bins_equal(result, expected), 'enforce_minimum_size does not match expected');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 6: enforce_minimum_size satisfies thresholds
% ---------------------------------------------------------------------------
fprintf('Test 6: enforce_minimum_size satisfies thresholds ... ');
obs    = load_observations(fixtures_dir);
pooled = mapa_calibrate(obs);
result = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
if height(result) > 1
    assert(all(result.n_obs  >= MIN_OBS),  'n_obs threshold violated');
    assert(all(result.n_bads >= MIN_BADS), 'n_bads threshold violated');
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 7: enforce_minimum_size preserves totals and monotonicity
% ---------------------------------------------------------------------------
fprintf('Test 7: enforce_minimum_size preserves totals and monotonicity ... ');
obs    = load_observations(fixtures_dir);
pooled = mapa_calibrate(obs);
result = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
assert(sum(result.n_obs)  == sum(pooled.n_obs),  'n_obs total mismatch after min size');
assert(sum(result.n_bads) == sum(pooled.n_bads), 'n_bads total mismatch after min size');
rates = result.n_bads ./ result.n_obs;
assert(all(diff(rates) <= 0), 'Bad rates not non-increasing after min size');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 8: enforce_minimum_size is noop with default thresholds
% ---------------------------------------------------------------------------
fprintf('Test 8: enforce_minimum_size is noop with default thresholds ... ');
obs    = load_observations(fixtures_dir);
pooled = mapa_calibrate(obs);
result = enforce_minimum_size(pooled);
assert(bins_equal(result, pooled), 'enforce_minimum_size with defaults changed bins');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 9: Bayesian adjustment matches expected
% ---------------------------------------------------------------------------
fprintf('Test 9: Bayesian adjustment matches expected ... ');
obs      = load_observations(fixtures_dir);
pooled   = mapa_calibrate(obs);
sized    = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
result   = apply_bayesian_adjustment(sized, BAYESIAN_K);
expected = load_bins(fixtures_dir, 'expected_calibrated_bins.csv');
expected.pd = double(expected.pd);
assert(height(result) == height(expected), 'Row count mismatch in calibrated bins');
assert(all(result.score_min == expected.score_min), 'score_min mismatch');
assert(all(result.score_max == expected.score_max), 'score_max mismatch');
assert(all(result.n_obs     == expected.n_obs),     'n_obs mismatch');
assert(all(result.n_bads    == expected.n_bads),    'n_bads mismatch');
for ii = 1:height(result)
    assert(abs(result.pd(ii) - expected.pd(ii)) < 1e-9, ...
        sprintf('pd mismatch at row %d: got %.15g, expected %.15g', ii, result.pd(ii), expected.pd(ii)));
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 10: Bayesian adjustment shrinks toward prior
% ---------------------------------------------------------------------------
fprintf('Test 10: Bayesian adjustment shrinks toward prior ... ');
obs    = load_observations(fixtures_dir);
pooled = mapa_calibrate(obs);
prior  = sum(pooled.n_bads) / sum(pooled.n_obs);
result = apply_bayesian_adjustment(pooled, BAYESIAN_K, prior);
for ii = 1:height(pooled)
    orig_rate = pooled.n_bads(ii) / pooled.n_obs(ii);
    adj_pd    = result.pd(ii);
    lo_bound  = min(orig_rate, prior);
    hi_bound  = max(orig_rate, prior);
    eps_tol   = 1e-12;
    assert(adj_pd >= lo_bound - eps_tol && adj_pd <= hi_bound + eps_tol, ...
        sprintf('pd at row %d not between bad_rate and prior', ii));
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 11: repool_calibrated_bins matches expected
% ---------------------------------------------------------------------------
fprintf('Test 11: repool_calibrated_bins matches expected ... ');
obs        = load_observations(fixtures_dir);
pooled     = mapa_calibrate(obs);
sized      = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
calibrated = apply_bayesian_adjustment(sized, BAYESIAN_K);
result     = repool_calibrated_bins(calibrated);
expected   = load_bins(fixtures_dir, 'expected_repooled_calibrated_bins.csv');
expected.pd = double(expected.pd);
assert(height(result) == height(expected), 'Row count mismatch in repooled bins');
assert(all(result.score_min == expected.score_min), 'score_min mismatch');
assert(all(result.score_max == expected.score_max), 'score_max mismatch');
assert(all(result.n_obs     == expected.n_obs),     'n_obs mismatch');
assert(all(result.n_bads    == expected.n_bads),    'n_bads mismatch');
for ii = 1:height(result)
    assert(abs(result.pd(ii) - expected.pd(ii)) < 1e-9, ...
        sprintf('pd mismatch at row %d: got %.15g, expected %.15g', ii, result.pd(ii), expected.pd(ii)));
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 12: repool_calibrated_bins restores monotonicity
% ---------------------------------------------------------------------------
fprintf('Test 12: repool_calibrated_bins restores monotonicity ... ');
obs        = load_observations(fixtures_dir);
pooled     = mapa_calibrate(obs);
sized      = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
calibrated = apply_bayesian_adjustment(sized, BAYESIAN_K);
% The fixture deliberately has a non-monotone pd sequence after Bayesian adjustment
pds_before = calibrated.pd;
assert(~all(diff(pds_before) <= 0), 'Expected non-monotone pds before repooling');
result     = repool_calibrated_bins(calibrated);
pds_after  = result.pd;
assert(all(diff(pds_after) <= 0), 'pds not non-increasing after repooling');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 13: repool_calibrated_bins preserves totals
% ---------------------------------------------------------------------------
fprintf('Test 13: repool_calibrated_bins preserves totals ... ');
obs        = load_observations(fixtures_dir);
pooled     = mapa_calibrate(obs);
sized      = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
calibrated = apply_bayesian_adjustment(sized, BAYESIAN_K);
result     = repool_calibrated_bins(calibrated);
assert(sum(result.n_obs)  == sum(calibrated.n_obs),  'n_obs total mismatch after repool');
assert(sum(result.n_bads) == sum(calibrated.n_bads), 'n_bads total mismatch after repool');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 14: interpolate_pd matches expected smoothed PDs
% ---------------------------------------------------------------------------
fprintf('Test 14: interpolate_pd matches expected smoothed PDs ... ');
obs        = load_observations(fixtures_dir);
pooled     = mapa_calibrate(obs);
sized      = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
calibrated = apply_bayesian_adjustment(sized, BAYESIAN_K);
repooled   = repool_calibrated_bins(calibrated);
smooth_exp = read_csv_table(fullfile(fixtures_dir, 'expected_smoothed_pds.csv'));
for ii = 1:height(smooth_exp)
    score  = double(smooth_exp.score(ii));
    exp_pd = double(smooth_exp.pd(ii));
    res    = interpolate_pd(repooled, score);
    rel_err = abs(res - exp_pd) / max(abs(exp_pd), 1e-15);
    assert(abs(res - exp_pd) < 1e-9 || rel_err < 1e-9, ...
        sprintf('interpolate_pd mismatch at score %g: got %.15g, expected %.15g', score, res, exp_pd));
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 15: interpolate_pd is monotone non-increasing
% ---------------------------------------------------------------------------
fprintf('Test 15: interpolate_pd is monotone non-increasing ... ');
obs        = load_observations(fixtures_dir);
pooled     = mapa_calibrate(obs);
sized      = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
calibrated = apply_bayesian_adjustment(sized, BAYESIAN_K);
repooled   = repool_calibrated_bins(calibrated);
scores     = unique([repooled.score_min; repooled.score_max], 'sorted');
pds        = arrayfun(@(s) interpolate_pd(repooled, s), scores);
assert(all(diff(pds) <= 1e-15), 'interpolate_pd not monotone non-increasing');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 16: interpolate_pd at pool midpoint matches pool pd
% ---------------------------------------------------------------------------
fprintf('Test 16: interpolate_pd at pool midpoint matches pool pd ... ');
obs        = load_observations(fixtures_dir);
pooled     = mapa_calibrate(obs);
sized      = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
calibrated = apply_bayesian_adjustment(sized, BAYESIAN_K);
repooled   = repool_calibrated_bins(calibrated);
for ii = 1:height(repooled)
    midpoint = (repooled.score_min(ii) + repooled.score_max(ii)) / 2;
    res      = interpolate_pd(repooled, midpoint);
    assert(abs(res - repooled.pd(ii)) < 1e-9, ...
        sprintf('interpolate_pd at midpoint mismatch for bin %d', ii));
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 17: run_pipeline bands match repool_calibrated_bins
% ---------------------------------------------------------------------------
fprintf('Test 17: run_pipeline bands match repool_calibrated_bins ... ');
obs        = load_observations(fixtures_dir);
pooled     = mapa_calibrate(obs);
sized      = enforce_minimum_size(pooled, MIN_OBS, MIN_BADS);
calibrated = apply_bayesian_adjustment(sized, BAYESIAN_K);
repooled   = repool_calibrated_bins(calibrated);
pipeline   = run_pipeline(obs, BAYESIAN_K, MIN_OBS, MIN_BADS);
assert(calibrated_bins_equal(pipeline.bands, repooled), ...
    'run_pipeline bands do not match repool_calibrated_bins output');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 18: run_pipeline pd_for_score matches interpolate_pd and expected
% ---------------------------------------------------------------------------
fprintf('Test 18: run_pipeline pd_for_score matches interpolate_pd and expected ... ');
obs      = load_observations(fixtures_dir);
pipeline = run_pipeline(obs, BAYESIAN_K, MIN_OBS, MIN_BADS);
smooth_exp2 = read_csv_table(fullfile(fixtures_dir, 'expected_smoothed_pds.csv'));
for ii = 1:height(smooth_exp2)
    score  = double(smooth_exp2.score(ii));
    exp_pd = double(smooth_exp2.pd(ii));
    from_fn     = pipeline.pd_for_score(score);
    from_interp = interpolate_pd(pipeline.bands, score);
    assert(abs(from_fn - from_interp) < 1e-15, ...
        sprintf('pd_for_score != interpolate_pd at score %g', score));
    rel_err = abs(from_fn - exp_pd) / max(abs(exp_pd), 1e-15);
    assert(abs(from_fn - exp_pd) < 1e-9 || rel_err < 1e-9, ...
        sprintf('pd_for_score mismatch at score %g', score));
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 19: mapa_pool min_confidence matches expected
% ---------------------------------------------------------------------------
fprintf('Test 19: mapa_pool min_confidence matches expected ... ');
obs      = load_observations(fixtures_dir);
initial  = bins_from_observations(obs);
result   = mapa_pool(initial, false, MIN_CONFIDENCE);
expected = load_bins(fixtures_dir, 'expected_pooled_bins_confidence.csv');
assert(bins_equal(result, expected), 'mapa_pool min_confidence does not match expected');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 20: mapa_pool min_confidence preserves totals and monotonicity
% ---------------------------------------------------------------------------
fprintf('Test 20: mapa_pool min_confidence preserves totals and monotonicity ... ');
obs     = load_observations(fixtures_dir);
initial = bins_from_observations(obs);
result  = mapa_pool(initial, false, MIN_CONFIDENCE);
assert(sum(result.n_obs)  == sum(initial.n_obs),  'n_obs total mismatch (confidence)');
assert(sum(result.n_bads) == sum(initial.n_bads), 'n_bads total mismatch (confidence)');
rates = result.n_bads ./ result.n_obs;
assert(all(diff(rates) <= 0), 'Bad rates not non-increasing (confidence)');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 21: mapa_pool min_confidence merges at least as much as plain mapa
% ---------------------------------------------------------------------------
fprintf('Test 21: mapa_pool confidence merges at least as much as plain ... ');
obs       = load_observations(fixtures_dir);
initial   = bins_from_observations(obs);
plain     = mapa_pool(initial);
confident = mapa_pool(initial, false, MIN_CONFIDENCE);
assert(height(confident) <= height(plain), ...
    'Confidence-based mapa did not merge at least as many bins as plain mapa');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 22: increasing direction
% ---------------------------------------------------------------------------
fprintf('Test 22: increasing direction ... ');
obs     = load_observations(fixtures_dir);
flipped = table(-obs.score, obs.bad, 'VariableNames', {'score', 'bad'});
result  = mapa_calibrate(flipped, true);
rates   = result.n_bads ./ result.n_obs;
assert(all(diff(rates) >= 0), 'Bad rates not non-decreasing in increasing mode');
assert(sum(result.n_obs) == height(obs), 'n_obs total mismatch in increasing mode');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 23: weighted bins_from_observations matches expected
% ---------------------------------------------------------------------------
fprintf('Test 23: weighted bins_from_observations matches expected ... ');
w_obs    = load_weighted_observations(fixtures_dir);
w_result = bins_from_observations(w_obs);
w_exp    = load_weighted_bins(fixtures_dir, 'expected_initial_bins_weighted.csv');
assert(weighted_bins_equal(w_result, w_exp), ...
    'weighted bins_from_observations does not match expected');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 24: weighted bins differ from unweighted
% ---------------------------------------------------------------------------
fprintf('Test 24: weighted bins differ from unweighted ... ');
w_obs       = load_weighted_observations(fixtures_dir);
w_bins      = bins_from_observations(w_obs);
uw_obs      = table(w_obs.score, w_obs.bad, 'VariableNames', {'score', 'bad'});
uw_bins     = bins_from_observations(uw_obs);
w_rates     = w_bins.n_bads ./ w_bins.n_obs;
uw_rates    = uw_bins.n_bads ./ uw_bins.n_obs;
common_len  = min(height(w_bins), height(uw_bins));
assert(any(abs(w_rates(1:common_len) - uw_rates(1:common_len)) > 1e-12), ...
    'Weighted and unweighted bad rates should differ');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 25: weighted pooling preserves totals
% ---------------------------------------------------------------------------
fprintf('Test 25: weighted pooling preserves totals ... ');
w_obs     = load_weighted_observations(fixtures_dir);
w_initial = bins_from_observations(w_obs);
w_pooled  = mapa_pool(w_initial);
assert(abs(sum(w_pooled.n_obs)  - sum(w_initial.n_obs))  < 1e-6, ...
    'weighted n_obs total mismatch after pooling');
assert(abs(sum(w_pooled.n_bads) - sum(w_initial.n_bads)) < 1e-6, ...
    'weighted n_bads total mismatch after pooling');
assert(sum(w_pooled.count)      == sum(w_initial.count), ...
    'weighted count total mismatch after pooling');
assert(sum(w_pooled.count_bads) == sum(w_initial.count_bads), ...
    'weighted count_bads total mismatch after pooling');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 26: weighted pooling is monotone
% ---------------------------------------------------------------------------
fprintf('Test 26: weighted pooling is monotone ... ');
w_obs    = load_weighted_observations(fixtures_dir);
w_pooled = mapa_pool(bins_from_observations(w_obs));
w_rates  = w_pooled.n_bads ./ w_pooled.n_obs;
assert(all(diff(w_rates) <= 0), ...
    'Weighted pooled bad rates are not non-increasing');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 27: weighted pooling matches expected
% ---------------------------------------------------------------------------
fprintf('Test 27: weighted pooling matches expected ... ');
w_obs    = load_weighted_observations(fixtures_dir);
w_pooled = mapa_pool(bins_from_observations(w_obs));
w_exp    = load_weighted_bins(fixtures_dir, 'expected_pooled_bins_weighted.csv');
assert(weighted_bins_equal(w_pooled, w_exp), ...
    'weighted pooling does not match expected');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 28: weighted enforce_minimum_size uses counts
% ---------------------------------------------------------------------------
fprintf('Test 28: weighted enforce_minimum_size uses counts ... ');
w_obs    = load_weighted_observations(fixtures_dir);
w_pooled = mapa_pool(bins_from_observations(w_obs));
w_sized  = enforce_minimum_size(w_pooled, MIN_OBS, MIN_BADS, false, [], true);
if height(w_sized) > 1
    assert(all(w_sized.count      >= MIN_OBS),  'count threshold violated (weighted)');
    assert(all(w_sized.count_bads >= MIN_BADS), 'count_bads threshold violated (weighted)');
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 29: weighted enforce_minimum_size matches expected
% ---------------------------------------------------------------------------
fprintf('Test 29: weighted enforce_minimum_size matches expected ... ');
w_obs    = load_weighted_observations(fixtures_dir);
w_pooled = mapa_pool(bins_from_observations(w_obs));
w_sized  = enforce_minimum_size(w_pooled, MIN_OBS, MIN_BADS, false, [], true);
w_exp    = load_weighted_bins(fixtures_dir, 'expected_min_size_bins_weighted.csv');
assert(weighted_bins_equal(w_sized, w_exp), ...
    'weighted enforce_minimum_size does not match expected');
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 30: weighted run_pipeline matches expected
% ---------------------------------------------------------------------------
fprintf('Test 30: weighted run_pipeline matches expected ... ');
w_obs      = load_weighted_observations(fixtures_dir);
w_pipeline = run_pipeline(w_obs, BAYESIAN_K, MIN_OBS, MIN_BADS, [], false, [], true);
w_exp      = load_weighted_bins(fixtures_dir, 'expected_repooled_calibrated_bins_weighted.csv');
assert(height(w_pipeline.bands) == height(w_exp), ...
    'Row count mismatch in weighted pipeline bands');
assert(all(w_pipeline.bands.score_min == w_exp.score_min), ...
    'score_min mismatch in weighted pipeline');
assert(all(w_pipeline.bands.score_max == w_exp.score_max), ...
    'score_max mismatch in weighted pipeline');
assert(all(abs(w_pipeline.bands.n_obs  - w_exp.n_obs)  < 1e-9), ...
    'n_obs mismatch in weighted pipeline');
assert(all(abs(w_pipeline.bands.n_bads - w_exp.n_bads) < 1e-9), ...
    'n_bads mismatch in weighted pipeline');
assert(all(w_pipeline.bands.count      == w_exp.count), ...
    'count mismatch in weighted pipeline');
assert(all(w_pipeline.bands.count_bads == w_exp.count_bads), ...
    'count_bads mismatch in weighted pipeline');
for ii = 1:height(w_exp)
    assert(abs(w_pipeline.bands.pd(ii) - w_exp.pd(ii)) < 1e-9, ...
        sprintf('pd mismatch at row %d in weighted pipeline', ii));
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
% Test 31: weighted smoothed PDs match expected
% ---------------------------------------------------------------------------
fprintf('Test 31: weighted smoothed PDs match expected ... ');
w_obs      = load_weighted_observations(fixtures_dir);
w_pipeline = run_pipeline(w_obs, BAYESIAN_K, MIN_OBS, MIN_BADS, [], false, [], true);
w_smooth   = read_csv_table(fullfile(fixtures_dir, 'expected_smoothed_pds_weighted.csv'));
for ii = 1:height(w_smooth)
    score  = double(w_smooth.score(ii));
    exp_pd = double(w_smooth.pd(ii));
    res_pd = w_pipeline.pd_for_score(score);
    rel_err = abs(res_pd - exp_pd) / max(abs(exp_pd), 1e-15);
    assert(abs(res_pd - exp_pd) < 1e-9 || rel_err < 1e-9, ...
        sprintf('weighted smoothed pd mismatch at score %g: got %.15g, expected %.15g', ...
                score, res_pd, exp_pd));
end
fprintf('PASS\n');

% ---------------------------------------------------------------------------
fprintf('\nAll tests passed.\n');

% Helper functions are in the private/ subdirectory for Octave compatibility.
