function k = linardKernDiagCompute(kern, x)

% LINARDKERNDIAGCOMPUTE Compute diagonal of linear ARD kernel.

% KERN

scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

k = sum(x.*x, 2)*kern.variance;
