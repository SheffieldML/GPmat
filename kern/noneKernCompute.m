function k = noneKernCompute(kern, x, x2)


% NONEKERNCOMPUTE Compute the NONE kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the dummy kernel function
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the dummy kernel function
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : noneKernParamInit, kernCompute, kernCreate, noneKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2008

% KERN



if nargin < 3
  k = spalloc(size(x, 1), size(x, 1), 0);
else
  k = spalloc(size(x, 1), size(x2, 1), 0);
end
