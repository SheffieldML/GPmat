function k = mlpKernDiagCompute(kern, x)

% MLPKERNDIAGCOMPUTE Compute diagonal of multi-layer perceptron kernel.

% KERN

% KERN


numer = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance;
denom = numer+1;
k = kern.variance*asin(numer./denom);
