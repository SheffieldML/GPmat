function k = biasKernDiagCompute(x, kern)

% BIASKERNDIAGCOMPUTE Compute diagonal of bias kernel.

% IVM

k = repmat(kern.variance, size(x, 1), 1);