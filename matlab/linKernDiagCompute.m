function k = linKernDiagCompute(kern, x)

% LINKERNDIAGCOMPUTE Compute diagonal of linear kernel.

% KERN

% KERN


linPart = ones(size(x, 1), 1);
k =  sum(x.*x, 2)*kern.variance;
