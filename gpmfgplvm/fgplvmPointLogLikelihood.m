function ll = fgplvmPointLogLikelihood(model, x, y)

% FGPLVMPOINTLOGLIKELIHOOD Log-likelihood of a point for the GP-LVM.
% FORMAT
% DESC returns the log likelihood of a latent point and an observed
% data point for the posterior prediction of the GP-LVM model.
% ARG model : the model for which the point prediction will be
% made.
% ARG x : the latent point for which the posterior distribution
% will be evaluated.
% ARG y : the observed data point for which the posterior
% distribution will be evaluated.
%
% SEEALSO : gpPointLogLikelihood, fgplvmCreate, fgplvmOptimisePoint, fgplvmPointObjective
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% FGPLVM

ll = gpPointLogLikelihood(model, x, y);
% check if there is a prior over latent space 
if isfield(model, 'prior')
  for i = 1:size(x, 1)
    ll = ll + priorLogProb(model.prior, x(i, :));
  end
end
