function g = modelLogLikeGradients(model)

% MODELLOGLIKEGRADIENTS Compute a model's gradients wrt log likelihood.
% FORMAT
% DESC is a wrapper function to compute the gradients of the log
% likelihood of a given model.
% ARG model : the model for which likelihoods are computed.
% RETURN g : teh gradients of the likelihood with respect to the
% parameters.
%
% SEEALSO : modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2005

% MLTOOLS

fhandle = str2func([model.type 'LogLikeGradients']);
g = fhandle(model);

if isfield(model, 'paramGroups')
  g = g*model.paramGroups;
end
