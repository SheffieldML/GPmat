function gX = cmpndKernDiagGradX(kern, x)

% CMPNDKERNDIAGGRADX Gradient of compound kernel's diagonal with respect to a point x.

% KERN


gX = kernDiagGradX(kern.comp{1}, x);
for i = 2:length(kern.comp)
  gX = gX + kernDiagGradX(kern.comp{i}, x);
end
