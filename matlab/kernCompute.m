function k = kernCompute(x, kern, x2)

% KERNCOMPUTE Compute the kernel given the parameters and X.

% IVM

if nargin < 3
  k = feval([kern.type 'KernCompute'], x, kern);
else
  k = feval([kern.type 'KernCompute'], x, kern, x2);
end
