function k = cmpndKernCompute(x, kern, x2)

% CMPNDKERNCOMPUTE Compute the kernel given the parameters and X.

% IVM

if nargin > 2
  k = kernCompute(x, kern.comp{1}, x2);
  for i = 2:length(kern.comp)
    k  = k + kernCompute(x, kern.comp{i}, x2);
  end
else
  k  = kernCompute(x, kern.comp{1});
  for i = 2:length(kern.comp)
    k  = k + kernCompute(x, kern.comp{i});
  end
end
