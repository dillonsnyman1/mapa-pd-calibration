function result = run_pipeline(observations, k, min_obs, min_bads, prior, increasing, min_confidence)
% RUN_PIPELINE  Run the full MAPA pipeline.
%
%   result = RUN_PIPELINE(observations, k)
%   result = RUN_PIPELINE(observations, k, min_obs, min_bads)
%   result = RUN_PIPELINE(observations, k, min_obs, min_bads, prior)
%   result = RUN_PIPELINE(observations, k, min_obs, min_bads, prior, increasing)
%   result = RUN_PIPELINE(observations, k, min_obs, min_bads, prior, increasing, min_confidence)
%
%   Chains:
%     bins_from_observations -> mapa_pool -> enforce_minimum_size ->
%     apply_bayesian_adjustment -> repool_calibrated_bins
%
%   observations   — table with `score` and `bad` columns, or N-by-2 matrix.
%   k              — Bayesian credibility weight; see apply_bayesian_adjustment.
%   min_obs        — minimum observations per bin (default 0).
%   min_bads       — minimum bads per bin (default 0).
%   prior          — PD to shrink toward (default [] = overall bad rate).
%   increasing     — logical (default false); direction of monotonicity.
%   min_confidence — optional confidence level for z-test merging.
%
%   Returns a struct with fields:
%     .bands       — table of final calibrated bins
%     .pd_for_score — function handle @(score) interpolate_pd(bands, score)

if nargin < 3 || isempty(min_obs);         min_obs        = 0;     end
if nargin < 4 || isempty(min_bads);        min_bads       = 0;     end
if nargin < 5;                             prior          = [];    end
if nargin < 6 || isempty(increasing);      increasing     = false; end
if nargin < 7;                             min_confidence = [];    end

pooled     = mapa_calibrate(observations, increasing, min_confidence);
sized      = enforce_minimum_size(pooled, min_obs, min_bads, increasing, min_confidence);
calibrated = apply_bayesian_adjustment(sized, k, prior);
bands      = repool_calibrated_bins(calibrated, increasing);

result.bands        = bands;
result.pd_for_score = @(score) interpolate_pd(bands, score);
end
