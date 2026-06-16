function bins = bins_from_observations(observations)
% BINS_FROM_OBSERVATIONS  Group raw (score, bad) observations into one bin
%                         per unique score, ordered by score ascending.
%
%   bins = BINS_FROM_OBSERVATIONS(observations)
%
%   observations — table with variables `score` (double) and `bad` (double,
%                  1 for a default and 0 otherwise), or an N-by-2 numeric
%                  matrix [score, bad].
%
%   bins — table with variables score_min, score_max, n_obs, n_bads,
%           one row per unique score value.

if isnumeric(observations)
    scores = observations(:, 1);
    bads   = observations(:, 2);
else
    scores = observations.score;
    bads   = observations.bad;
end

unique_scores = unique(scores, 'sorted');
n = numel(unique_scores);

score_min_col = zeros(n, 1);
score_max_col = zeros(n, 1);
n_obs_col     = zeros(n, 1);
n_bads_col    = zeros(n, 1);

for i = 1:n
    s = unique_scores(i);
    idx = scores == s;
    score_min_col(i) = s;
    score_max_col(i) = s;
    n_obs_col(i)     = sum(idx);
    n_bads_col(i)    = sum(bads(idx));
end

bins = table(score_min_col, score_max_col, n_obs_col, n_bads_col, ...
    'VariableNames', {'score_min', 'score_max', 'n_obs', 'n_bads'});
end
