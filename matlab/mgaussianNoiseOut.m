function y = mgaussianNoiseOut(noise, mu, varsigma)

% MGAUSSIANNOISEOUT Ouput from Variable variance Gaussian noise model.

% NOISE

D = size(mu, 2);
y = zeros(size(mu));
for i = 1:D
  y(:, i) = mu(:, i) + noise.bias(i);
end

