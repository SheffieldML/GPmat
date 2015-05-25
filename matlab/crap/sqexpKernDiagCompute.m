function k = sqexpKernDiagCompute(x, kern)

% SQEXPKERNDIAGCOMPUTE Compute diagonal of squared exponential kernel.

% IVM

rbfPart = ones(size(x, 1), 1);
k = rbfPart*(kern.rbfVariance + kern.whiteVariance) + kern.biasVariance;
