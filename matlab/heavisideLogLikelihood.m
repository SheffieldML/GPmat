function L = heavisideLogLikelihood(noise, mu, varsigma, y)

% HEAVISIDELOGLIKELIHOOD Log-likelihood of data under heaviside noise model.

% IVM

L = sum(sum(log(heavisideLikelihood(noise, mu, varsigma, y))));
