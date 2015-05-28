function [k, rbfPart, n2] = sqexpKernCompute(kern, x, x2)


% SQEXPKERNCOMPUTE Compute the SQEXP kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the pre-built compound squared exponential
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the pre-built compound squared exponential
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : sqexpKernParamInit, kernCompute, kernCreate, sqexpKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004

% KERN


if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  k = rbfPart + kern.whiteVariance*eye(size(x, 1));
else
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  k = rbfPart;
end
k = k + kern.biasVariance;
