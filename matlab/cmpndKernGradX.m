function gX = cmpndKernGradX(kern, x, X2)

% CMPNDKERNGRADX Gradient of compound kernel with respect to a point x.

% KERN


gX = kernGradX(kern.comp{1}, x, X2);
for i = 2:length(kern.comp)
  gX = gX + kernGradX(kern.comp{i}, x, X2);
end
