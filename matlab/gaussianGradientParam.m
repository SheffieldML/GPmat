function g = gaussianGradientParam(model, params)

% GAUSSIANGRADIENTPARAM Gradient of the Gaussian noise model's parameters.

% IVM

D = size(model.y, 2);
u = zeros(size(model.y));

for i = 1:D
  mu(:, i) = model.mu(:, i) + model.noise.bias(i);
end

u = model.y - mu;
nu = 1./(model.varSigma+model.noise.sigma2);
u = u.*nu;
gbias = sum(u, 1);
gsigma2 = sum(nu - u.*u, 1);
gsigma2 = -0.5*gsigma2.*model.noise.sigma2;
g = [gbias sum(gsigma2)];