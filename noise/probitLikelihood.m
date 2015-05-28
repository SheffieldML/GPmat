function L = probitLikelihood(noise, mu, varsigma, y)

% PROBITLIKELIHOOD Likelihood of data under probit noise model.

% NOISE

D = size(y, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
L = cumGaussian((y.*mu)./(sqrt(noise.sigma2+varsigma)));
