function bins = bins_from_observations(observations)
% BINS_FROM_OBSERVATIONS  Group raw (score, bad) observations into one bin
%                         per unique score, ordered by score ascending.
%
%   bins = BINS_FROM_OBSERVATIONS(observations)
%
%   observations — table with variables `score` (double) and `bad` (double,
%                  1 for a default and 0 otherwise), or an N-by-2 numeric
%                  matrix [score, bad].  An optional third column (or table
%                  variable `weight`) supplies per-observation weights.
%                  When omitted, weight = 1 for every row.
%
%   bins — table with variables score_min, score_max, n_obs, n_bads,
%          count, count_bads — one row per unique score value.
%          n_obs / n_bads are weighted sums; count / count_bads are raw
%          observation counts.  When no weight column is supplied the two
%          pairs are identical.

if isnumeric(observations)
    scores = observations(:, 1);
    bads   = observations(:, 2);
    if size(observations, 2) >= 3
        weights = observations(:, 3);
    else
        weights = ones(size(scores));
    end
else
    scores = observations.score;
    bads   = observations.bad;
    if any(strcmp(observations.Properties.VariableNames, 'weight'))
        weights = observations.weight;
    else
        weights = ones(size(scores));
    end
end

unique_scores = unique(scores, 'sorted');
n = numel(unique_scores);

score_min_col  = zeros(n, 1);
score_max_col  = zeros(n, 1);
n_obs_col      = zeros(n, 1);
n_bads_col     = zeros(n, 1);
count_col      = zeros(n, 1);
count_bads_col = zeros(n, 1);

for i = 1:n
    s = unique_scores(i);
    idx = scores == s;
    score_min_col(i)  = s;
    score_max_col(i)  = s;
    n_obs_col(i)      = sum(weights(idx));
    n_bads_col(i)     = sum(bads(idx) .* weights(idx));
    count_col(i)      = sum(idx);
    count_bads_col(i) = sum(bads(idx));
end

bins = table(score_min_col, score_max_col, n_obs_col, n_bads_col, ...
    count_col, count_bads_col, ...
    'VariableNames', {'score_min', 'score_max', 'n_obs', 'n_bads', ...
                      'count', 'count_bads'});
end
