function k = rbfKernDiagCompute(kern, x)

% RBFKERNDIAGCOMPUTE Compute diagonal of rbf kernel.

% KERN

% KERN


rbfPart = ones(size(x, 1), 1);
k = rbfPart*kern.variance;
