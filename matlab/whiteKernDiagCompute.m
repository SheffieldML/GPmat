function k = whiteKernDiagCompute(kern, x)

% WHITEKERNDIAGCOMPUTE Compute diagonal of white noise kernel.

% KERN


k = repmat(kern.variance, size(x, 1), 1);