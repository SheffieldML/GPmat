function g = modelOutputGrad(model, X, dim)

% MODELOUTPUTGRAD Compute derivatives with respect to params of model outputs.
% FORMAT
% DESC gives the gradients of the outputs from the model with
% respect to the parameters for a given set of inputs.
% ARG model : the model structure for which gradients are computed.
% ARG X : input locations where gradients are to be computed.
% RETURN g : gradients of the model output with respect to the
% model parameters for the given input locations. The size of the
% returned matrix is of dimension number of data x number of
% parameters x number of model outputs (which maintains
% compatability with NETLAB).
%
% FORMAT
% DESC gives the gradients of the outputs from the model with
% respect to the parameters for a given set of inputs.
% ARG model : the model structure for which gradients are computed.
% ARG X : input locations where gradients are to be comxfputed.
% ARG dim : the dimension of the model for which gradients are required.
% RETURN g : gradients of the model output with respect to the
% model parameters for the given input locations. The size of the
% returned matrix is of dimension number of data x number of
% parameters.
%
% SEEALSO : modelCreate, modelLogLikelihood, modelLogLikeGradients, mlpderiv
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Cark Henrik Ek, 2007

% MLTOOLS

if(nargin>2)
  fhandle = str2func([model.type 'OutputGrad']);
  g = fhandle(model, X, dim);
else
  fhandle = str2func([model.type 'OutputGrad']);
  gtemp = fhandle(model, X);
  
  if isfield(model, 'paramGroups')
    g = zeros(size(X, 1), size(model.paramGroups, 2), size(gtemp, 3));
    for i = 1:size(gtemp, 3)
      g = gtemp(:, :, i)*model.paramGroups;
    end
  else 
    g = gtemp;
  end
end
