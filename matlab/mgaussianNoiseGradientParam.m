function g = mgaussianNoiseGradientParam(noise, mu, varsigma, y)

% MGAUSSIANNOISEGRADIENTPARAM Gradient of the Variable variance Gaussian noise model's parameters.

% IVM


D = size(y, 2);
u = zeros(size(y));
nu = zeros(size(y));
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
  nu(:, i) = 1./(varsigma(:, i) + noise.sigma2(i));
end

u = y - mu;

% Remove unlabelled points from the gradient.
u(find(isnan(y))) = 0;
nu(find(isnan(y)))= 0;
u = u.*nu;
gbias = sum(u, 1);
gsigma2 = -.5*sum(nu - u.*u, 1);
g = [gbias gsigma2];
