function g = modelOutputGrad(model, X)

% MODELOUTPUTGRAD Compute derivatives with respect to params of model outputs.
% FORMAT
% DESC gives the gradients of the outputs from the model with
% respect to the parameters for a given set of inputs.
% ARG model : the model structure for which gradients are computed.
% ARG X : input locations where gradients are to be computed.
% RETURN g : gradients of the model output with respect to the
% model parameters for the given input locations.
%
% SEEALSO : modelCreate, modelLogLikelihood, modelLogLikeGradients
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS

fhandle = str2func([model.type 'OutputGrad']);
gtemp = fhandle(model, X);

if isfield(model, 'paramGroups')
  g = zeros(size(X, 1), size(kern.paramGroups, 2), size(gtemp, 3));
  for i = 1:size(gtemp, 3)
    g = gtemp(:, :, i)*kern.paramGroups;
  end
else 
  g = gtemp;
end