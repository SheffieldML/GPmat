function k = sqexpKernDiagCompute(kern, x)

% SQEXPKERNDIAGCOMPUTE Compute diagonal of squared exponential kernel.

% KERN


rbfPart = ones(size(x, 1), 1);
k = rbfPart*(kern.rbfVariance + kern.whiteVariance) + kern.biasVariance;
