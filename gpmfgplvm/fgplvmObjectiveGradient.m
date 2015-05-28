function [f, g] = fgplvmObjectiveGradient(params, model)

% FGPLVMOBJECTIVEGRADIENT Wrapper function for FGPLVM objective and gradient.
% FORMAT
% DESC returns the negative log likelihood of a Gaussian process
% model given the model structure and a vector of parameters. This
% allows the use of NETLAB minimisation functions to find the model
% parameters.
% ARG params : the parameters of the model for which the objective
% will be evaluated.
% ARG model : the model structure for which the objective will be
% evaluated.
% RETURN f : the negative log likelihood of the FGPLVM model.
% RETURN g : the gradient of the negative log likelihood of the FGPLVM
% model with respect to the parameters.
%
% SEEALSO : minimize, fgplvmCreate, fgplvmGradient, fgplvmLogLikelihood, fgplvmOptimise
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% FGPLVM

% Check how the optimiser has given the parameters
if size(params, 1) > size(params, 2)
  % As a column vector ... transpose everything.
  transpose = true;
  model = fgplvmExpandParam(model, params');
else
  transpose = false;
  model = fgplvmExpandParam(model, params);
end

f = - fgplvmLogLikelihood(model);
if nargout > 1
  g = - fgplvmLogLikeGradients(model);
end
if transpose
  g = g';
end
