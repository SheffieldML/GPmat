function L = orderedLogLikelihood(noise, mu, varsigma, y)

% ORDEREDLOGLIKELIHOOD Log-likelihood of data under ordered categorical noise model.

% IVM

L = sum(sum(log(orderedLikelihood(noise, mu, varsigma, y))));