function k = rbfardKernDiagCompute(kern, x)

% RBFARDKERNDIAGCOMPUTE Compute diagonal of radial basis function ARD kernel.

% KERN

k = repmat(kern.variance, size(x, 1), 1);
