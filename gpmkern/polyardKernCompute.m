function [k, innerProd, arg, denom, numer, vec] = polyardKernCompute(kern, x, x2)


% POLYARDKERNCOMPUTE Compute the POLYARD kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the automatic relevance determination polynomial
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the automatic relevance determination polynomial
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : polyardKernParamInit, kernCompute, kernCreate, polyardKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

if nargin < 3
  innerProd = x*x';
  arg = innerProd*kern.weightVariance + kern.biasVariance;
  k = kern.variance*arg.^kern.degree;
else
  x2 = x2*scales;
  innerProd = x*x2';  
  arg = innerProd*kern.weightVariance + kern.biasVariance;
  k = kern.variance*arg.^kern.degree;
end
