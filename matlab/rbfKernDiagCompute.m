function k = rbfKernDiagCompute(kern, x)

% RBFKERNDIAGCOMPUTE Compute diagonal of rbf kernel.

% KERN

k = repmat(kern.variance, size(x, 1), 1);
