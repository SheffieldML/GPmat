function k = rbfard2KernDiagCompute(kern, x)

% RBFARD2KERNDIAGCOMPUTE Compute diagonal of RBFARD2 kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the automatic relevance determination radial basis function kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : rbfard2KernParamInit, kernDiagCompute, kernCreate, rbfard2KernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% COPYRIGHT : Michalis K. Titsias, 2009

% KERN


k = repmat(kern.variance, size(x, 1), 1);
