function k = mlpardKernDiagCompute(kern, x)

% MLPARDKERNDIAGCOMPUTE Compute diagonal of multi-layer perceptron ARD kernel.

% KERN

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
numer = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance;
denom = numer+1;
k = kern.variance*asin(numer./denom);
