function k = linKernCompute(x, kern, x2)

% LINKERNCOMPUTE Compute the kernel given the parameters and X.

% IVM

if nargin < 3
   k = x*x'*kern.variance;
else
  k = x*x2'*kern.variance;
end
