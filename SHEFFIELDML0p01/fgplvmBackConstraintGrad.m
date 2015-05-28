function g = fgplvmBackConstraintGrad(model, gX)

% FGPLVMBACKCONSTRAINTGRAD Gradient with respect to back constraints if present.
%
%	Description:
%
%	FGPLVMBACKCONSTRAINTGRAD(MODEL, GX) converts the gradients of the
%	GP-LVM model log likelihood with respect to the latent positions to
%	be gradients with respect to the parameters of the back constraints.
%	 Arguments:
%	  MODEL - the GP-LVM model structure for which the conversion is to
%	   be done.
%	  GX - the gradients of the log likelihood with respect to the back
%	   constraint parameters.
%	
%
%	See also
%	% SEEALSO FGPLVMLOGLIKEGRADIENTS, FGPLVMCREATE


%	Copyright (c) 2005, 2006, 2007 Neil D. Lawrence


% Check for back constraints.
if isfield(model, 'back') & ~isempty(model.back)
  g_w = modelOutputGrad(model.back, model.y);
  g_modelParams = zeros(size(g_w, 2), 1);
  for i = 1:model.q
    g_modelParams = g_modelParams + g_w(:, :, i)'*gX(:, i);
  end
  g = g_modelParams;
else
  % Do nothing
  g = gX;
end