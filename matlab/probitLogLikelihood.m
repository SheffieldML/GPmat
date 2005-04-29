function L = probitLogLikelihood(noise, mu, varsigma, y)

% PROBITLOGLIKELIHOOD Log-likelihood of data under probit noise model.

% NOISE

D = size(y, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end

L = sum(sum(lnCumGaussian((y.*mu)./(sqrt(noise.sigma2+varsigma)))));
