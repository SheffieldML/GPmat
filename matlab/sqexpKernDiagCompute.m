function k = sqexpKernDiagCompute(kern, x)

% SQEXPKERNDIAGCOMPUTE Compute diagonal of squared exponential kernel.

% KERN

k = repmat(kern.rbfVariance+kern.whiteVariance, size(x, 1), 1)  + kern.biasVariance;
