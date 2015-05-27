function ll = modelLogLikelihood(model)

% MODELLOGLIKELIHOOD Compute a model log likelihood.
% FORMAT
% DESC computes the log likelihood of the given model.
% ARG model : the model for which the log likelihood is to be
% computed.
% RETURN ll : the log likelihood of the data given the model.
%
% SEEALSO : modelLogLikeGradients, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS

fhandle = str2func([model.type 'LogLikelihood']);
ll = fhandle(model);
