function ll = fgplvmLogLikelihood(model)

% FGPLVMLOGLIKELIHOOD Log-likelihood for a GP-LVM.

% FGPLVM


ll = gpLogLikelihood(model);

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  % A dynamics model is being used.
  ll = ll + modelLogLikelihood(model.dynamics);
elseif isfield(model, 'prior') &  ~isempty(model.prior)
  for i = 1:model.N
    ll = ll + priorLogProb(model.prior, model.X(i, :));
  end
end

switch model.approx
  case {'dtc', 'fitc', 'pitc'}
   if isfield(model, 'inducingPrior') &  ~isempty(model.inducingPrior)
     for i = 1:model.k
       ll = ll + priorLogProb(model.inducingPrior, model.X_u(i, :));    
     end
   end
 otherwise
  % do nothing
end


