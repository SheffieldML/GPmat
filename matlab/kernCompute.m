function k = kernCompute(kern, x, x2)

% KERNCOMPUTE Compute the kernel given the parameters and X.

% IVM

if nargin < 3
  k = feval([kern.type 'KernCompute'], kern, x);
else
  k = feval([kern.type 'KernCompute'], kern, x, x2);
end
