function ll = fgplvmPointLogLikelihood(model, x, y)

% FGPLVMPOINTLOGLIKELIHOOD Log-likelihood of a point for the GP-LVM.
%
% ll = fgplvmPointLogLikelihood(model, x, y)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmPointLogLikelihood.m version 1.2



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
