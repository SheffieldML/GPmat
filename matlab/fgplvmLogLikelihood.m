function ll = fgplvmLogLikelihood(model)

% FGPLVMLOGLIKELIHOOD Log-likelihood for a GP-LVM.

% FGPLVM


ll = gpLogLikelihood(model);

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  % Dynamics kernel is being used.
  ll = ll + gpLogLikelihood(model.dynamics, model.X(2:end, :));
  % Use prior on first point only.
  if isfield(model, 'prior') &  ~isempty(model.prior)
    ll = ll + priorLogProb(model.prior, model.X(1, :));
  end  
elseif isfield(model, 'prior') &  ~isempty(model.prior)
  for i = 1:model.N
    ll = ll + priorLogProb(model.prior, model.X(i, :));
  end
end


