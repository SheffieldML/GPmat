function gX = cmpndKernDiagGradX(x, kern)

% CMPNDKERNDIAGGRADX Gradient of compound kernel's diagonal with respect to a point x.

% IVM

gX = kernDiagGradX(x, kern.comp{1});
for i = 2:length(kern.comp)
  gX = gX + kernDiagGradX(x, kern.comp{i});
end
