function ll = fgplvmLogLikelihood(model)

% FGPLVMLOGLIKELIHOOD Log-likelihood for a GP-LVM.
% FORMAT
% DESC returns the log likelihood for a given GP-LVM model.
% ARG model : the model for which the log likelihood is to be
% computed. The model contains the data for which the likelihood is
% being computed in the 'y' component of the structure.
% RETURN ll : the log likelihood of the data given the model.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009
%
% SEEALSO : gpLogLikelihood, fgplvmCreate
%
% MODIFICATIONS : Carl Henrik Ek, 2008, 2009
%
% FGPLVM

 
ll = gpLogLikelihood(model);

if isfield(model, 'dynamics') && ~isempty(model.dynamics)
  % A dynamics model is being used.
  ll = ll + modelLogLikelihood(model.dynamics);
elseif isfield(model, 'prior') &&  ~isempty(model.prior)
  for i = 1:model.N
    ll = ll + priorLogProb(model.prior, model.X(i, :));
  end
end

switch model.approx
  case {'dtc', 'dtcvar', 'fitc', 'pitc'}
   if isfield(model, 'inducingPrior') && ~isempty(model.inducingPrior)
     for i = 1:model.k
       ll = ll + priorLogProb(model.inducingPrior, model.X_u(i, :));    
     end
   end
 otherwise
  % do nothing
end

if(isfield(model,'constraints')&&~isempty(model.constraints))
  for(i = 1:1:model.constraints.numConstraints)
    ll = ll + constraintLogLikelihood(model.constraints.comp{i},model.X);
  end
end

