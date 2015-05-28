function [m, beta] = scaleNoiseSites(noise, g, nu, mu, varSigma, y)

% SCALENOISESITES Site updates for Scale noise model.

% NOISE

N = size(y, 1);
D = length(noise.bias);
beta = zeros(N, D);
for i = 1:size(y, 2)
  m(:, i) = (y(:, i) - noise.bias(i))/noise.scale(i);
end
beta = repmat(1./noise.sigma2, N, D);
