function y = gaussianNoiseOut(noise, mu, varsigma)

% GAUSSIANNOISEOUT Output from Gaussian noise model.

% IVM

D = size(mu, 2);
y = zeros(size(mu));
for i = 1:D
  y(:, i) = mu(:, i) + noise.bias(i);
end
