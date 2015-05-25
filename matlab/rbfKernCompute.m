function [k, sk, n2] = rbfKernCompute(kern, x, x2)

% RBFKERNCOMPUTE Compute the RBF kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the radial basis function
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only the exponential part).
%
% FORMAT
% DESC computes the kernel matrix for the radial basis function
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
% RETURN sk : unscaled kernel matrix (i.e. only the exponential part).
%
% SEEALSO : rbfKernParamInit, kernCompute, kernCreate, rbfKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% MODIFICATIONS : Mauricio Alvarez, 2009, David Luengo, 2009

% KERN

if nargin < 3
  n2 = dist2(x, x);
else
  n2 = dist2(x, x2);
end

wi2 = (.5 .* kern.inverseWidth);
sk = exp(-n2*wi2);
k = kern.variance*sk;
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    k = k * sqrt(kern.inverseWidth/(2*pi));
end
