function k = indexKernDiagCompute(kern, x)

% INDEXKERNDIAGCOMPUTE Compute diagonal of INDEX kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the index based covariance function kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : indexKernParamInit, kernDiagCompute, kernCreate, indexKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN
  k = repmat(kern.variance, size(x, 1), 1);
end
