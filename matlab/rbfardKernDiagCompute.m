function k = rbfardKernDiagCompute(kern, x)

% RBFARDKERNDIAGCOMPUTE Compute diagonal of radial basis function ARD kernel.

% IVM


k = ones(size(x, 1), 1)*kern.variance;
