function gX = cmpndKernGradX(x, kern, X2)

% CMPNDKERNGRADX Gradient of compound kernel with respect to a point x.

% IVM

gX = kernGradX(x, kern.comp{1}, X2);
for i = 2:length(kern.comp)
  gX = gX + kernGradX(x, kern.comp{i}, X2);
end
