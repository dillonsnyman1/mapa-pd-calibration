function result = pd_violates(lower, upper, increasing)
% PD_VIOLATES  Whether upper (higher-scoring bin) violates monotonicity of
%              pd relative to lower.
%
%   result = PD_VIOLATES(lower, upper, increasing)
%
%   lower, upper  — single-row calibrated bin tables (with pd column)
%   increasing    — logical; if true, pd must be non-decreasing

if increasing
    result = upper.pd < lower.pd;
else
    result = upper.pd > lower.pd;
end
end
