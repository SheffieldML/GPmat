function L = probitLogLikelihood(noise, mu, varsigma, y)

% PROBITLOGLIKELIHOOD Log-likelihood of data under probit noise model.

% IVM

L = sum(sum(log(probitLikelihood(noise, mu, varsigma, y))));
