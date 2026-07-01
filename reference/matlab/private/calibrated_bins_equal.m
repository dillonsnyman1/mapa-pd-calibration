function ok = calibrated_bins_equal(a, b, tol)
    if nargin < 3; tol = 1e-9; end
    if ~bins_equal(a, b); ok = false; return; end
    rel_err = abs(a.pd - b.pd) ./ max(abs(b.pd), 1e-15);
    ok = all(abs(a.pd - b.pd) < tol | rel_err < tol);
end
