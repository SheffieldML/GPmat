function ll = fgplvmPointLogLikelihood(model, x, y)

% FGPLVMPOINTLOGLIKELIHOOD Log-likelihood of a point for the GP-LVM.

% FGPLVM

logTwoPi = log(2*pi);
[mu, varSigma] = gpPosteriorMeanVar(model, x);
ll = 0;
ydiff = y-mu;
ll = log(varSigma) + (ydiff.*ydiff)./varSigma +logTwoPi;
ll(find(isnan(ll)))=0;
ll = -0.5*sum(ll, 2);
% check if there is a prior over latent space 
if isfield(model, 'prior')
  ll = ll + priorLogProb(model.prior, x);
end
