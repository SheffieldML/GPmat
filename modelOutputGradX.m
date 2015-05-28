function g = modelOutputGradX(model, X)

% MODELOUTPUTGRADX Compute derivatives with respect to model inputs of model outputs.
% FORMAT
% DESC gives the gradients of the outputs from the model with
% respect to the inputs.
% ARG model : the model structure for which gradients are computed.
% ARG X : input locations where gradients are to be computed.
% RETURN g : gradients of the model output with respect to the
% model parameters for the given input locations.
%
% SEEALSO : modelCreate, modelOutputGrad, modelLogLikelihood, modelLogLikeGradients
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

fhandle = str2func([model.type 'OutputGradX']);
g = fhandle(model, X);

