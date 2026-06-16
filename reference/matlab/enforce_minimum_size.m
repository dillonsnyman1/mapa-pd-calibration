function bins = enforce_minimum_size(bins, min_obs, min_bads, increasing, min_confidence)
% ENFORCE_MINIMUM_SIZE  Pool bins below minimum size thresholds, then re-run
%                       mapa_pool to restore monotonicity.
%
%   bins = ENFORCE_MINIMUM_SIZE(bins)
%   bins = ENFORCE_MINIMUM_SIZE(bins, min_obs, min_bads)
%   bins = ENFORCE_MINIMUM_SIZE(bins, min_obs, min_bads, increasing)
%   bins = ENFORCE_MINIMUM_SIZE(bins, min_obs, min_bads, increasing, min_confidence)
%
%   bins           — table of bins (score_min, score_max, n_obs, n_bads),
%                    typically from mapa_pool.
%   min_obs        — minimum observations per bin (default 0).
%   min_bads       — minimum bads per bin (default 0).
%   increasing     — logical (default false); forwarded to mapa_pool.
%   min_confidence — optional confidence level; forwarded to mapa_pool.
%
%   A bin violates if n_obs < min_obs or n_bads < min_bads. Each violating
%   bin is merged into the adjacent bin with the closer bad rate.

if nargin < 2 || isempty(min_obs);  min_obs  = 0; end
if nargin < 3 || isempty(min_bads); min_bads = 0; end
if nargin < 4 || isempty(increasing);       increasing     = false; end
if nargin < 5;                              min_confidence = [];    end

% Work row-by-row; convert to cell array of single-row tables for easy splicing
n = height(bins);
bin_list = cell(n, 1);
for i = 1:n
    bin_list{i} = bins(i, :);
end

while numel(bin_list) > 1
    n = numel(bin_list);

    % Find first violating bin
    violator = 0;
    for i = 1:n
        b = bin_list{i};
        if b.n_obs < min_obs || b.n_bads < min_bads
            violator = i;
            break;
        end
    end

    if violator == 0
        break;
    end

    if violator == 1
        neighbour = 2;
    elseif violator == n
        neighbour = n - 1;
    else
        rate       = bin_list{violator}.n_bads / bin_list{violator}.n_obs;
        left_diff  = abs(rate - bin_list{violator-1}.n_bads / bin_list{violator-1}.n_obs);
        right_diff = abs(rate - bin_list{violator+1}.n_bads / bin_list{violator+1}.n_obs);
        if left_diff <= right_diff
            neighbour = violator - 1;
        else
            neighbour = violator + 1;
        end
    end

    lo = min(violator, neighbour);
    hi = max(violator, neighbour);

    merged = merge_bins(bin_list{lo}, bin_list{hi});

    % Rebuild list: before lo, merged, after hi
    before = bin_list(1:lo-1);
    after  = bin_list(hi+1:end);
    bin_list = [before; {merged}; after];
end

bins = vertcat(bin_list{:});
bins = mapa_pool(bins, increasing, min_confidence);
end
