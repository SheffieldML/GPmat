function g = ivmGradX(model, x, y);

% IVMGRADX Returns the gradient of the log-likelihood wrt x.
%
%	Description:
%
%	IVMGRADX(MODEL, X, Y) returns the gradient of the approximate log
%	likelihood with respect to an input location x. This is used for
%	optimising with respect to x in the GP-LVM.
%	 Arguments:
%	  MODEL - the model for which the gradient is being computed.
%	  X - the input location where the gradient is to be evaluated.
%	  Y - the target location where the gradient is being evaluated.
%	
%
%	See also
%	% SEEALSO IVMPOSTERIORMEANVAR, IVMPOSTERIORGRADMEANVAR


%	Copyright (c) 2004, 2005 Neil D. Lawrence


[mu, varsigma] = ivmPosteriorMeanVar(model, x);
[dmu, dvs] = ivmPosteriorGradMeanVar(model, x);

g = noiseGradX(model.noise, mu, varsigma, dmu, dvs, y);
