function k = whiteKernDiagCompute(x, kern)

% WHITEKERNDIAGCOMPUTE Compute diagonal of white noise kernel.

% IVM

k = repmat(kern.variance, size(x, 1), 1);