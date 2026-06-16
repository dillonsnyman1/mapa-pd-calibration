function bins = mapa_pool(bins, increasing, min_confidence)
% MAPA_POOL  Run the Monotone Adjacent Pooling Algorithm (PAVA-style).
%
%   bins = MAPA_POOL(bins)
%   bins = MAPA_POOL(bins, increasing)
%   bins = MAPA_POOL(bins, increasing, min_confidence)
%
%   bins           — table with variables score_min, score_max, n_obs,
%                    n_bads; ordered by score ascending (e.g. from
%                    bins_from_observations).
%   increasing     — logical (default false). If false, bad rate must be
%                    non-increasing. If true, non-decreasing.
%   min_confidence — optional scalar in (0,1). Adjacent bins whose bad rates
%                    are not distinguishable at this confidence level
%                    (two-proportion z-test) are merged even if they don't
%                    violate monotonicity. Omit or pass [] to disable.
%
%   Returns pooled bins as a table.

if nargin < 2 || isempty(increasing)
    increasing = false;
end
if nargin < 3
    min_confidence = [];
end

use_confidence = ~isempty(min_confidence);

% Stack: cell array of single-row tables
stack = {};

for i = 1:height(bins)
    b = bins(i, :);
    stack{end+1} = b; %#ok<AGROW>

    while numel(stack) >= 2
        lower = stack{end-1};
        upper = stack{end};

        should_merge = bin_violates(lower, upper, increasing) || ...
            (use_confidence && not_significant(lower, upper, min_confidence));

        if ~should_merge
            break;
        end

        merged = merge_bins(lower, upper);
        stack(end-1:end) = [];
        stack{end+1} = merged; %#ok<AGROW>
    end
end

bins = vertcat(stack{:});
end
