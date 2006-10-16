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
% SEEALSO : fgplvmCreate, fgplvmOptimisePoint, fgplvmPointObjective
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% FGPLVM

logTwoPi = log(2*pi);
[mu, varSigma] = gpPosteriorMeanVar(model, x);
ll = zeros(size(x, 1), 1);
ydiff = y-mu;
ll = log(varSigma) + (ydiff.*ydiff)./varSigma +logTwoPi;
ll(find(isnan(ll)))=0;
ll = -0.5*sum(ll, 2);
% check if there is a prior over latent space 
if isfield(model, 'prior')
  for i = 1:size(x, 1)
    ll = ll + priorLogProb(model.prior, x(i, :));
  end
end
