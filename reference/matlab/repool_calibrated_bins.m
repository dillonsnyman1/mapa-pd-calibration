function bins = repool_calibrated_bins(bins, increasing)
% REPOOL_CALIBRATED_BINS  Re-pool Bayesian-adjusted bins to restore
%                         monotonicity of pd.
%
%   bins = REPOOL_CALIBRATED_BINS(bins)
%   bins = REPOOL_CALIBRATED_BINS(bins, increasing)
%
%   bins       — table with columns score_min, score_max, n_obs, n_bads, pd;
%                typically from apply_bayesian_adjustment.
%   increasing — logical (default false). If false, pd must be non-increasing.
%                If true, non-decreasing.
%
%   Merging uses n_obs-weighted average of pd.

if nargin < 2 || isempty(increasing)
    increasing = false;
end

% Stack: cell array of single-row calibrated bin tables
stack = {};

for i = 1:height(bins)
    b = bins(i, :);
    stack{end+1} = b; %#ok<AGROW>

    while numel(stack) >= 2
        lower = stack{end-1};
        upper = stack{end};

        if ~pd_violates(lower, upper, increasing)
            break;
        end

        merged = merge_calibrated(lower, upper);
        stack(end-1:end) = [];
        stack{end+1} = merged; %#ok<AGROW>
    end
end

bins = vertcat(stack{:});
end
