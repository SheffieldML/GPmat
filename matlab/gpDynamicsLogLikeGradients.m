function g = gpDynamicsLogLikeGradients(model)

% GPDYNAMICSLOGLIKEGRADIENTS Gradients of the GP dynamics wrt parameters.
% FORMAT
% DESC Computes the gradients with respect to the log likelihood of
% the GP dynamics in a GP-LVM model.
% ARG model : the GP model for which log likelihood is to be
% computed.
% RETURN g : the gradients of the log likelihood with respect to
% the latent points and (optionally) parameters.
%
% SEEALSO : gpLogLikeGradients, gpDynamicsCreate, gpDynamicsLogLikelihood, modelLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009
%
% MODIFICATIONS : Carl Henrik Ek, 2006
  
% FGPLVM

if model.k ==0 & ~model.learn & ~model.learnScales
  g = [];
  return
end

g = gpLogLikeGradients(model);

if ~model.learn
  % If we aren't learning model parameters extract only X_u;
  % this is inefficient (but neater in the code) as we have also computed parameters 
  if ~model.learnScales
    if isfield(model, 'fixInducing') & model.fixInducing
      g = [];
    else
      g = g(1:model.k*model.q);
    end
  else
    switch model.approx
     case 'ftc'
      g =  [g(end-model.d + 1:end)];
     case {'dtc', 'dtcvar', 'fitc', 'pitc'}
      if isfield(model, 'fixInducing') & model.fixInducing
        g =  g(end-model.d:end-1);
      else
        g =  [g(1:model.k*model.q) g(end-model.d:end-1)];
      end
    end
  end
end
