function ll = gpReversibleDynamicsLogLikelihood(model)

% GPREVERSIBLEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the dynamics part.

% FGPLVM

ll = gpLogLikelihood(model);

% Use prior on first two points only for reversible dynamics model.
if isfield(model, 'prior') &  ~isempty(model.prior)
  ll = ll + priorLogProb(model.prior, model.X(1, :));
  ll = ll + priorLogProb(model.prior, model.X(2, :));
end  
