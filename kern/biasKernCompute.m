function [k, sk] = biasKernCompute(kern, x, x2)


% BIASKERNCOMPUTE Compute the BIAS kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the bias
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the inpute matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the bias
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : biasKernParamInit, kernCompute, kernCreate, biasKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


if nargin< 3
  dims = [size(x, 1), size(x, 1)];
else
  dims = [size(x, 1), size(x2, 1)];
end
k = repmat(kern.variance, dims);
if nargout > 1
  sk = ones(dims);
end
