function k = kernCompute(kern, x, x2)

% KERNCOMPUTE Compute the kernel given the parameters and X.

% KERN

fhandle = str2func([kern.type 'KernCompute']);
if nargin < 3
  k = fhandle(kern, x);
else
  k = fhandle(kern, x, x2);
end
