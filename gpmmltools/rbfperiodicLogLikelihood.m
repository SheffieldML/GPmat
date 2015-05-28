function ll = rbfperiodicLogLikelihood(model)

% RBFPERIODICLOGLIKELIHOOD Log likelihood of RBFPERIODIC model.
% FORMAT
% DESC computes the log likelihood of  the periodic radial basis function model.
% ARG model : the model structure for which log likelihood is being computed.
% RETURN ll : the computed log likelihood.
%
% SEEALSO : rbfperiodicCreate, rbfperiodicLogLikeGradients, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

ypred = rbfperiodicOut(model, model.X);

ll = - 0.5.*(sum(sum((model.y-ypred).^2)));
ll = real(ll);
