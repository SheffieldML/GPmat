function [k, rbfPart, n2] = sqexpKernCompute(x, kern, x2)

% SQEXPKERNCOMPUTE Compute the squared exponential kernel given the parameters and X.

% IVM


if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  k = rbfPart + kern.whiteVariance*eye(size(x, 1));
else
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  rbfPart = kern.rbfVariance*exp(-n2*wi2);
  k = rbfPart;
end
k = k + kern.biasVariance;