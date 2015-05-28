function k = polyardKernDiagCompute(kern, x)


% POLYARDKERNDIAGCOMPUTE Compute diagonal of POLYARD kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the automatic relevance determination polynomial kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : polyardKernParamInit, kernDiagCompute, kernCreate, polyardKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
arg = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance;
k = kern.variance*arg.^kern.degree;
