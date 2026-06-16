function bins = mapa_calibrate(observations, increasing, min_confidence)
% MAPA_CALIBRATE  Convenience wrapper: bins_from_observations + mapa_pool.
%
%   bins = MAPA_CALIBRATE(observations)
%   bins = MAPA_CALIBRATE(observations, increasing)
%   bins = MAPA_CALIBRATE(observations, increasing, min_confidence)
%
%   observations   — table with variables `score` and `bad`, or N-by-2
%                    numeric matrix [score, bad].
%   increasing     — logical (default false); see mapa_pool.
%   min_confidence — optional confidence level; see mapa_pool.
%
%   Returns pooled bins as a table.

if nargin < 2 || isempty(increasing)
    increasing = false;
end
if nargin < 3
    min_confidence = [];
end

bins = mapa_pool(bins_from_observations(observations), increasing, min_confidence);
end
