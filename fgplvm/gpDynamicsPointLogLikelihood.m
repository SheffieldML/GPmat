function ll = gpDynamicsPointLogLikelihood(model, x, y)

% GPDYNAMICSPOINTLOGLIKELIHOOD Compute the log likelihood of a point under the GP dynamics model.
% DESC computes the log likelihood of a given point under the GP
% dynamics model.
% ARG model : the model for which the log likelihood is to be
% computed.
% ARG x : the input location for the model.
% ARG y : the target location for the model.
% RETURN ll : the log likelihood of the given point.
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% SEEALSO : gpDynamicsCreate, gpDynamicsLogLikelihood, gpPointLogLikelihood
% 

% MODIFICATION: Carl Henrik Ek, 2009

% FGPLVM

if(isfield(model,'indexIn'))
  x = x(:,model.indexIn);
  y = y(:,model.indexIn);
end
ll = gpPointLogLikelihood(model, x, y);
