function [k, sk, n2] = ratquadKernCompute(kern, x, x2)

% RATQUADKERNCOMPUTE Compute the RATQUAD kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the rational quadratic
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the rational quadratic
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : ratquadKernParamInit, kernCompute, kernCreate, ratquadKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009

% KERN

wi2 = .5/(kern.lengthScale*kern.lengthScale*kern.alpha);
if nargin < 3
  n2 = dist2(x, x);
else
  n2 = dist2(x, x2);
end
sk = (1+n2*wi2).^-kern.alpha;
k = kern.variance*sk;
