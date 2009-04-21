function k = mlpKernDiagCompute(kern, x)


% MLPKERNDIAGCOMPUTE Compute diagonal of MLP kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the multi-layer perceptron kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : mlpKernParamInit, kernDiagCompute, kernCreate, mlpKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


numer = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance;
denom = numer+1;
k = 2/pi*kern.variance*asin(numer./denom);
