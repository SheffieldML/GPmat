function g = probitNoiseGradientParam(noise, mu, varsigma, y)

% PROBITNOISEGRADIENTPARAM Gradient of the probit noise model's parameters.

% NOISE

c = y./sqrt(noise.sigma2 + varsigma);
for i = 1:size(mu, 2)
  u(:, i) = c(:, i).*(mu(:, i) + noise.bias(i));
end
g = sum(c.*gradLogCumGaussian(u), 1);