function L = probitLogLikelihood(X, Y, model)

% PROBITLOGLIKELIHOOD Log-likelihood of data under probit noise model.

% IVM

L = sum(sum(log(probitLikelihood(X, Y, model))));
