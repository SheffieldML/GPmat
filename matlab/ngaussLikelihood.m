function L = ngaussLikelihood(noise, mu, varsigma, y)

% NGAUSSLIKELIHOOD Likelihood of data under noiseless Gaussian noise model.

% IVM

L = gaussianLikelihood(noise, mu, varsigma, y);