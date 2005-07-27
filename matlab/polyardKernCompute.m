function [k, innerProd, arg, denom, numer, vec] = polyardKernCompute(kern, x, x2)

% POLYARDKERNCOMPUTE Compute the polynomial ARD kernel given the parameters and X.

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
