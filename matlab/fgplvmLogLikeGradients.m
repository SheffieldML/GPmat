function g = fgplvmLogLikeGradients(model)

% FGPLVMLOGLIKEGRADIENTS Compute the gradients for the FGPLVM.
% FORMAT
% DESC returns the gradients of the log likelihood with respect to the
% parameters of the GP-LVM model and with respect to the latent
% positions of the GP-LVM model. 
% ARG model : the FGPLVM structure containing the parameters and
% the latent positions.
% RETURN g : the gradients of the latent positions and the
% parameters of the GP-LVM model.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : fgplvmLogLikelihood, fgplvmCreate, modelLogLikeGradients

% FGPLVM


[gParam, gX_u, gX] = gpLogLikeGradients(model);

gDynParam = [];
% Check if Dynamics kernel is being used.
if isfield(model, 'dynamics') & ~isempty(model.dynamics)

  % Get the dynamics parameters
  gDynParam = modelLogLikeGradients(model.dynamics);
  
  % Include the dynamics latent gradients.
  gX = gX + modelLatentGradients(model.dynamics);

elseif isfield(model, 'prior') &  ~isempty(model.prior)
  gX = gX + priorGradient(model.prior, model.X); 
end

switch model.approx
 case {'dtc', 'fitc', 'pitc'}
  if isfield(model, 'inducingPrior') & ~isempty(model.inducingPrior)
    gX_u = gX_u + priorGradient(model.inducingPrior, model.X_u);
  end
 otherwise
  % do nothing
end

% Concatanate existing parameter gradients.
gParam = [gParam gDynParam];

% Decide where to include gX_u.
if ~strcmp(model.approx, 'ftc') & model.fixInducing
  gX(model.inducingIndices, :) = gX(model.inducingIndices, :) + gX_u;
else
  gParam = [gX_u(:)' gParam];
end

% Check for back constraints.
if isfield(model, 'back')
  g_w = modelOutputGrad(model.back, model.y);
  g_modelParams = zeros(size(g_w, 2), 1);
  for i = 1:model.q
    g_modelParams = g_modelParams + g_w(:, :, i)'*gX(:, i);
  end
  g = [g_modelParams(:)' gParam];
else
  g = [gX(:)' gParam];
end




