function [g, nu] = gaussianNoiseNuG(noise, mu, varSigma, y)

% GAUSSIANNOISENUG Update nu and g parameters associated with Gaussian noise model.

% NOISE

D = size(y, 2);
nu = 1./(noise.sigma2+varSigma);
g = zeros(size(nu));
for i = 1:D
  g(:, i) = y(:, i) - mu(:, i) - ...
      noise.bias(i);
end
g = g.*nu;
