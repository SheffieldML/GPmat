function [k, n2] = rbfKernCompute(x, kern, x2)

% RBFKERNCOMPUTE Compute the kernel given the parameters and X.

% IVM

if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
else
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
end
