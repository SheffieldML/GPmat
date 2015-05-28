function ll = fgplvmPointLogLikelihood(model, x, y)

% FGPLVMPOINTLOGLIKELIHOOD Log-likelihood of a point for the GP-LVM.
%
%	Description:
%
%	FGPLVMPOINTLOGLIKELIHOOD(MODEL, X, Y) returns the log likelihood of
%	a latent point and an observed data point for the posterior
%	prediction of the GP-LVM model.
%	 Arguments:
%	  MODEL - the model for which the point prediction will be made.
%	  X - the latent point for which the posterior distribution will be
%	   evaluated.
%	  Y - the observed data point for which the posterior distribution
%	   will be evaluated.
%	
%
%	See also
%	GPPOINTLOGLIKELIHOOD, FGPLVMCREATE, FGPLVMOPTIMISEPOINT, FGPLVMPOINTOBJECTIVE


%	Copyright (c) 2005, 2006 Neil D. Lawrence


ll = gpPointLogLikelihood(model, x, y);
% check if there is a prior over latent space 
if isfield(model, 'prior')
  for i = 1:size(x, 1)
    ll = ll + priorLogProb(model.prior, x(i, :));
  end
end
