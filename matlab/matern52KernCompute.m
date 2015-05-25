function [K, sk] = matern52KernCompute(kern, x, x2)

% MATERN52KERNCOMPUTE Compute the MATERN52 kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the matern kernel with nu=5/2
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the matern kernel with nu=5/2
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : matern52KernParamInit, kernCompute, kernCreate, matern52KernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009

% KERN

if nargin < 3
  n2 = dist2(x, x);
else
  n2 = dist2(x, x2);
end
wi2 = (5/(kern.lengthScale*kern.lengthScale));
n2wi2 = n2*wi2;
sqrtn2wi2 = sqrt(n2wi2);
sk = (1+sqrtn2wi2+n2wi2/3).*exp(-sqrtn2wi2);
K = kern.variance*sk;
