function L = heavisideLogLikelihood(X, Y, model)

% HEAVISIDELOGLIKELIHOOD Log-likelihood of data under heaviside noise model.

% IVM

L = sum(sum(log(heavisideLikelihood(X, Y, model))));
