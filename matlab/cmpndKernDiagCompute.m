function k = cmpndKernDiagCompute(kern, x)

% CMPNDKERNDIAGCOMPUTE Compute diagonal of compound kernel.

% KERN

% KERN

k = kernDiagCompute(kern.comp{1}, x);
for i = 2:length(kern.comp)
  k  = k + kernDiagCompute(kern.comp{i}, x);
end
