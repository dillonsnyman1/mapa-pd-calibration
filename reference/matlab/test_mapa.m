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
% Local helper functions are defined at the bottom of this file (required by
% MATLAB: script local functions must follow all executable statements).

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
        sprintf('pd mismatch at row %d', ii));
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
smooth_exp = readtable(fullfile(fixtures_dir, 'expected_smoothed_pds.csv'));
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
smooth_exp2 = readtable(fullfile(fixtures_dir, 'expected_smoothed_pds.csv'));
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
fprintf('\nAll tests passed.\n');

% ===========================================================================
% Local helper functions (must appear after all executable script statements)
% ===========================================================================

function bins = load_bins(fixtures_dir, filename)
% Load a bin table from a CSV fixture file.
    t = readtable(fullfile(fixtures_dir, filename));
    bins = table( ...
        double(t.score_min), ...
        double(t.score_max), ...
        double(t.n_obs), ...
        double(t.n_bads), ...
        'VariableNames', {'score_min', 'score_max', 'n_obs', 'n_bads'});
    if any(strcmp(t.Properties.VariableNames, 'pd'))
        bins.pd = double(t.pd);
    end
end

function obs = load_observations(fixtures_dir)
% Load raw observations from CSV.
    t = readtable(fullfile(fixtures_dir, 'raw_observations.csv'));
    obs = table(double(t.score), double(t.bad), 'VariableNames', {'score', 'bad'});
end

function ok = bins_equal(a, b)
% TRUE if two bin tables have the same shape, scores, and counts.
    ok = height(a) == height(b) && ...
         all(a.score_min == b.score_min) && ...
         all(a.score_max == b.score_max) && ...
         all(a.n_obs     == b.n_obs) && ...
         all(a.n_bads    == b.n_bads);
end

function ok = calibrated_bins_equal(a, b, tol)
% TRUE if two calibrated bin tables agree within tolerance on all columns.
    if nargin < 3; tol = 1e-9; end
    if ~bins_equal(a, b); ok = false; return; end
    rel_err = abs(a.pd - b.pd) ./ max(abs(b.pd), 1e-15);
    ok = all(abs(a.pd - b.pd) < tol | rel_err < tol);
end
