function L = orderedLogLikelihood(X, Y, model)

% ORDEREDLOGLIKELIHOOD Log-likelihood of data under Ordered categorical model.

% IVM

L = sum(sum(log(orderedLikelihood(X, Y, model))));