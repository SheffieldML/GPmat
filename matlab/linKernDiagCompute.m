function k = linKernDiagCompute(kern, x)

% LINKERNDIAGCOMPUTE Compute diagonal of linear kernel.

% KERN

k =  sum(x.*x, 2)*kern.variance;
