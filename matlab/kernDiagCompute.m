function k = kernDiagCompute(x, kern, x2)

% KERNELCOMPUTE Compute the kernel given the parameters and X.

% IVM

if nargin < 3
  k = feval([kern.type 'KernDiagCompute'], x, kern);
else
  k = feval([kern.type 'KernDiagCompute'], x, kern, x2);
end
