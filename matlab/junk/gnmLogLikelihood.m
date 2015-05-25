function L = gnmLogLikelihood(X, Y, model)

% GNMLOGLIKELIHOOD Log-likelihood of data under gap noise model.

% IVM

L = sum(sum(log(gnmLikelihood(X, Y, model))));