function k = sqexpKernDiagCompute(kern, x)


% SQEXPKERNDIAGCOMPUTE Compute diagonal of SQEXP kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the pre-built compound squared exponential kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : sqexpKernParamInit, kernDiagCompute, kernCreate, sqexpKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004

% KERN


k = repmat(kern.rbfVariance+kern.whiteVariance, size(x, 1), 1)  + kern.biasVariance;
