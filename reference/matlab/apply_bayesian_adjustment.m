function calibrated = apply_bayesian_adjustment(bins, k, prior)
% APPLY_BAYESIAN_ADJUSTMENT  Shrink each bin's bad rate toward a prior using
%                            Bayesian (credibility) weighting.
%
%   calibrated = APPLY_BAYESIAN_ADJUSTMENT(bins, k)
%   calibrated = APPLY_BAYESIAN_ADJUSTMENT(bins, k, prior)
%
%   bins  — table of bins (score_min, score_max, n_obs, n_bads), typically
%           from mapa_pool.
%   k     — credibility weight (equivalent observations of the prior).
%           Larger k means stronger shrinkage.
%   prior — PD to shrink toward (scalar). If omitted or [], defaults to the
%           overall bad rate: sum(n_bads) / sum(n_obs).
%
%   Each bin's adjusted PD:
%       pd = (n_bads + k * prior) / (n_obs + k)
%
%   Returns a table with columns score_min, score_max, n_obs, n_bads, pd.

if nargin < 3 || isempty(prior)
    prior = sum(bins.n_bads) / sum(bins.n_obs);
end

pd = (bins.n_bads + k * prior) ./ (bins.n_obs + k);

calibrated = [bins, table(pd)];
end
