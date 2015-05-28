function ll = fgplvmDynamicsLogLikelihood(model)

% FGPLVMDYNAMICSLOGLIKELIHOOD Gradients of the dynamics portion of the log likelihood..


if ~isfield(model, 'dynamics')
  error('This model does not contain dynamics');
end
model.dynamics.Y = model.X(2:end, :);
ll = ll + fgplvmLogLikelihood(model.dynamics);
  