function pd = interpolate_pd(bins, score)
% INTERPOLATE_PD  Return a smoothed PD for an individual score via log-odds
%                 interpolation between pool midpoints.
%
%   pd = INTERPOLATE_PD(bins, score)
%
%   bins  — table with columns score_min, score_max, n_obs, n_bads, pd;
%           typically from repool_calibrated_bins.
%   score — scalar score value.
%
%   Each pool is anchored at its midpoint (score_min + score_max) / 2.
%   log-odds is linearly interpolated between the two bracketing midpoints,
%   then converted back to a probability:
%       log_odds = log((1 - pd) / pd)
%       pd = 1 / (1 + exp(log_odds))
%
%   Flat extrapolation beyond the first and last midpoints.

n        = height(bins);
mids     = (bins.score_min + bins.score_max) / 2;
log_odds = log((1 - bins.pd) ./ bins.pd);

if score <= mids(1)
    pd = bins.pd(1);
    return;
end
if score >= mids(n)
    pd = bins.pd(n);
    return;
end

for i = 1:n-1
    if mids(i) <= score && score <= mids(i+1)
        t            = (score - mids(i)) / (mids(i+1) - mids(i));
        interpolated = log_odds(i) + t * (log_odds(i+1) - log_odds(i));
        pd           = 1 / (1 + exp(interpolated));
        return;
    end
end

error('interpolate_pd: unreachable — score must lie between mids(1) and mids(n)');
end
