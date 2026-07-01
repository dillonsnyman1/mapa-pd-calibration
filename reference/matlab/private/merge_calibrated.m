function merged = merge_calibrated(a, b)
% MERGE_CALIBRATED  Merge two adjacent calibrated bins using n_obs-weighted
%                   average of pd.
%
%   merged = MERGE_CALIBRATED(a, b)
%
%   a, b   — single-row calibrated bin tables (score_min, score_max,
%             n_obs, n_bads, pd, count, count_bads)
%   merged — single-row calibrated bin table

n_obs  = double(a.n_obs) + double(b.n_obs);
n_bads = double(a.n_bads) + double(b.n_bads);
pd     = (double(a.pd) * double(a.n_obs) + double(b.pd) * double(b.n_obs)) / n_obs;
merged = table( ...
    a.score_min, b.score_max, n_obs, n_bads, ...
    a.count + b.count, a.count_bads + b.count_bads, pd, ...
    'VariableNames', {'score_min', 'score_max', 'n_obs', 'n_bads', 'count', 'count_bads', 'pd'});
end
