function [k, innerProd, arg, denom, numer, vec] = mlpardKernCompute(kern, x, x2)

% MLPARDKERNCOMPUTE Compute the multi-layer perceptron ARD kernel given the parameters and X.

% KERN



scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;

if nargin < 3
  innerProd = x*x';
  numer = innerProd*kern.weightVariance + kern.biasVariance;
  vec = diag(numer) + 1;
  denom = sqrt(vec*vec');
  arg = numer./denom;
  k = kern.variance*asin(arg);
else
  x2 = x2*scales;
  innerProd = x*x2';  
  numer = innerProd*kern.weightVariance + kern.biasVariance;
  vec1 = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance + 1;
  vec2 = sum(x2.*x2, 2)*kern.weightVariance + kern.biasVariance + 1;
  denom = sqrt(vec1*vec2');
  arg = numer./denom;
  k = kern.variance*asin(arg);
end
