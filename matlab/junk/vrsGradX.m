function g = vrsGradX(X, Y, model)

% VRSGRADX Gradient wrt x of log-likelihood for various noise models.

% IVM

if size(X, 1) > 1
  error('This function only takes one data-point');
end


g = noiseGradX(x, y, model.noise.comp{1});
for i = 2:length(kern.comp)
  gX = gX + kernGradX(kern.comp{i}, x, X2);
end
