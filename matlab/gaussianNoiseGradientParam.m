function g = gaussianNoiseGradientParam(noise, mu, varsigma, y)

% GAUSSIANNOISEGRADIENTPARAM Gradient of the Gaussian noise model's parameters.

% IVM

D = size(y, 2);
u = zeros(size(y));

for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end

u = y - mu;
nu = 1./(varsigma+noise.sigma2);
u = u.*nu;
gbias = sum(u, 1);
gsigma2 = sum(nu - u.*u, 1);
gsigma2 = -0.5*gsigma2.*noise.sigma2;
g = [gbias sum(gsigma2)];