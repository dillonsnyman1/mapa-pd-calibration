function result = bin_violates(lower, upper, increasing)
% BIN_VIOLATES  Whether upper (higher-scoring bin) violates monotonicity of
%               bad_rate relative to lower.
%
%   result = BIN_VIOLATES(lower, upper, increasing)
%
%   lower, upper  — single-row bin tables (score_min, score_max, n_obs, n_bads)
%   increasing    — logical; if true, bad rate must be non-decreasing

lower_rate = lower.n_bads / lower.n_obs;
upper_rate = upper.n_bads / upper.n_obs;

if increasing
    result = upper_rate < lower_rate;
else
    result = upper_rate > lower_rate;
end
end
