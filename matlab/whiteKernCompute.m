function [k, sk] = whiteKernCompute(kern, x, x2)


% WHITEKERNCOMPUTE Compute the WHITE kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the white noise
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the inpute matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the white noise
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : whiteKernParamInit, kernCompute, kernCreate, whiteKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

% KERN

if nargin < 3
  sk = speye(size(x, 1));
  k = kern.variance*sk;
else
  k = spalloc(size(x, 1), size(x2, 1), 0);
end
