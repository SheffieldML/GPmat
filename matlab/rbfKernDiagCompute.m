function k = rbfKernDiagCompute(kern, x)

% RBFKERNDIAGCOMPUTE Compute diagonal of rbf kernel.

% IVM

rbfPart = ones(size(x, 1), 1);
k = rbfPart*kern.variance;
