function k = biasKernDiagCompute(kern, x)

% BIASKERNDIAGCOMPUTE Compute diagonal of bias kernel.

% KERN


k = repmat(kern.variance, size(x, 1), 1);