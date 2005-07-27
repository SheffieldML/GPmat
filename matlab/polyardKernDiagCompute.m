function k = polyardKernDiagCompute(kern, x)

% POLYARDKERNDIAGCOMPUTE Compute diagonal of multi-layer perceptron ARD kernel.

% KERN

scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
arg = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance;
k = kern.variance*arg.^kern.degree;
