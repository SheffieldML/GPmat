function [k, n2] = rbfardKernCompute(kern, x, x2)

% RBFARDKERNCOMPUTE Compute the radial basis function ARD kernel given the parameters and X.

% KERN



scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
    
if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
else
  x2 = x2*scales;
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
end
