function [nu, g] = gaussianNoiseUpdateParams(noise, mu, varsigma, y, index)

% GAUSSIANNOISEUPDATEPARAMS Update parameters for Gaussian noise model.

% IVM

D = size(y, 2);
nu = 1./(noise.sigma2+varsigma(index, :));
g = zeros(size(nu));
for i = 1:D
  g(:, i) = y(index, i) - mu(index, i) - ...
      noise.bias(i);
end
g = g.*nu;