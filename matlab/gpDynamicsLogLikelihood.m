function ll = gpDynamicsLogLikelihood(model)

% GPDYNAMICSLOGLIKELIHOOD Give the log likelihood of the dynamics part.

% FGPLVM

ll = gpLogLikelihood(model);

% Use prior on first point only for dynamics models.
if isfield(model, 'prior') &  ~isempty(model.prior)
  ll = ll + priorLogProb(model.prior, model.X(1, :));
end  
