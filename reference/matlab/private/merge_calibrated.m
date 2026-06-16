function merged = merge_calibrated(a, b)
% MERGE_CALIBRATED  Merge two adjacent calibrated bins using n_obs-weighted
%                   average of pd.
%
%   merged = MERGE_CALIBRATED(a, b)
%
%   a, b   — single-row calibrated bin tables (score_min, score_max,
%             n_obs, n_bads, pd)
%   merged — single-row calibrated bin table

n_obs  = a.n_obs + b.n_obs;
n_bads = a.n_bads + b.n_bads;
pd     = (a.pd * a.n_obs + b.pd * b.n_obs) / n_obs;

merged = table( ...
    a.score_min, ...
    b.score_max, ...
    n_obs, ...
    n_bads, ...
    pd, ...
    'VariableNames', {'score_min', 'score_max', 'n_obs', 'n_bads', 'pd'});
end
