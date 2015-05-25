function L = vrsLogLikelihood(X, Y, model)

% VRSLOGLIKELIHOOD Log-likelihood of data under various noise models.

% IVM

L = sum(sum(log(vrsLikelihood(X, Y, model))));