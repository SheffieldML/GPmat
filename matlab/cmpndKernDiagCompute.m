function k = cmpndKernDiagCompute(x, kern)

% CMPNDKERNDIAGCOMPUTE Compute diagonal of compound kernel.

% IVM
k = kernDiagCompute(x, kern.comp{1});
for i = 2:length(kern.comp)
  k  = k + kernDiagCompute(x, kern.comp{i});
end
