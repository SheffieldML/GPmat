function L = ngaussLogLikelihood(noise, mu, varsigma, y)

% NGAUSSLOGLIKELIHOOD Log-likelihood of data under noiseless Gaussian noise model.

% IVM

L = gaussianLogLikelihood(noise, mu, varsigma, y);
