function k = kernDiagCompute(kern, x, x2)

% KERNELCOMPUTE Compute the kernel given the parameters and X.

% KERN


if nargin < 3
  k = feval([kern.type 'KernDiagCompute'], kern, x);
else
  k = feval([kern.type 'KernDiagCompute'], kern, x, x2);
end
