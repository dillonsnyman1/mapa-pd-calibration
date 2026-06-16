function result = not_significant(a, b, confidence)
% NOT_SIGNIFICANT  Two-proportion z-test: returns true when the bad rates
%                  of bins a and b are NOT significantly different at the
%                  given confidence level (i.e. the two bins should be merged).
%
%   result = NOT_SIGNIFICANT(a, b, confidence)
%
%   a, b        — single-row bin tables (n_obs, n_bads)
%   confidence  — confidence level, e.g. 0.95
%
%   Uses sqrt(2)*erfinv(confidence) for the critical z-value, which is
%   equivalent to norminv((1+confidence)/2) but requires no Statistics
%   Toolbox — erfinv is a base MATLAB function.

pooled_rate = (a.n_bads + b.n_bads) / (a.n_obs + b.n_obs);

if pooled_rate <= 0 || pooled_rate >= 1
    result = true;
    return;
end

se         = sqrt(pooled_rate * (1 - pooled_rate) * (1/a.n_obs + 1/b.n_obs));
rate_a     = a.n_bads / a.n_obs;
rate_b     = b.n_bads / b.n_obs;
z          = abs(rate_a - rate_b) / se;
z_critical = sqrt(2) * erfinv(confidence);  % = norminv((1+confidence)/2)
result     = z < z_critical;
end
