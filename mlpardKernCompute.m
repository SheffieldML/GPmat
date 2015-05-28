function [k, sk, innerProd, arg, denom, numer, vec1, vec2] = mlpardKernCompute(kern, x, x2)


% MLPARDKERNCOMPUTE Compute the MLPARD kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the automatic relevance determination multi-layer perceptron
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the automatic relevance determination multi-layer perceptron
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : mlpardKernParamInit, kernCompute, kernCreate, mlpardKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

if nargin < 3
  innerProd = x*x';
  numer = innerProd*kern.weightVariance + kern.biasVariance;
  vec1 = diag(numer) + 1;
  denom = sqrt(vec1*vec1');
  arg = numer./denom;
  sk = 2/pi*asin(arg);
  k = kern.variance*sk;
else
  x2 = x2*scales;
  innerProd = x*x2';  
  numer = innerProd*kern.weightVariance + kern.biasVariance;
  vec1 = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance + 1;
  vec2 = sum(x2.*x2, 2)*kern.weightVariance + kern.biasVariance + 1;
  denom = sqrt(vec1*vec2');
  arg = numer./denom;
  sk = 2/pi*asin(arg);
  k = kern.variance*sk;
end
