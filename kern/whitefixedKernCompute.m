function [k, sk] = whitefixedKernCompute(kern, x, x2)


% WHITEFIXEDKERNCOMPUTE Compute the WHITEFIXED kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the fixed parameter white noise
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the inpute matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the fixed parameter white noise
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : whitefixedKernParamInit, kernCompute, kernCreate, whitefixedKernDiagCompute
%
% COPYRIGHT : Nathaniel J. King, 2006
%
% MODIFICATIONS : Neil D. Lawrence, 2009
  
% KERN


if nargin < 3
  k = whiteKernCompute(kern, x);
else
  k = whiteKernCompute(kern, x, x2);
end
if nargout > 1 
  sk = k;
end
