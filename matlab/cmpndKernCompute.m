function k = cmpndKernCompute(kern, x, x2)

% CMPNDKERNCOMPUTE Compute the kernel given the parameters and X.

% KERN

% KERN


if nargin > 2
  k = kernCompute(kern.comp{1}, x, x2);
  for i = 2:length(kern.comp)
    k  = k + kernCompute(kern.comp{i}, x, x2);
  end
else
  k  = kernCompute(kern.comp{1}, x);
  for i = 2:length(kern.comp)
    k  = k + kernCompute(kern.comp{i}, x);
  end
end
