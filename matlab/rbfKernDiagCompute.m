function k = rbfKernDiagCompute(x, kern)

% RBFKERNDIAGCOMPUTE Compute diagonal of rbf kernel.

% IVM

rbfPart = ones(size(x, 1), 1);
k = rbfPart*kern.variance;
