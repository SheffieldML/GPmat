function [m, beta] = gaussianNoiseSites(noise, g, nu, mu, varSigma, y)

% GAUSSIANNOISESITES Site updates for Gaussian noise model.

% IVM

N = size(y, 1);
D = length(noise.bias);
beta = zeros(N, D);
for i = 1:size(y, 2)
  m(:, i) = y(:, i) - noise.bias(i);
end
beta = repmat(1./noise.sigma2, N, D);