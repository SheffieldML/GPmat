function y = heavisideNoiseOut(noise, mu, varsigma)

% HEAVISIDENOISEOUT Output from heaviside noise model.

% IVM
D = size(mu, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
y = sign(mu);
