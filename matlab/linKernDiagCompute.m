function k = linKernDiagCompute(x, kern)

% LINKERNDIAGCOMPUTE Compute diagonal of linear kernel.

% IVM

linPart = ones(size(x, 1), 1);
k =  sum(x.*x, 2)*kern.variance;
