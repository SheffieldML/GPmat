function [mu, varsigma] = fgplvmPosteriorMeanVar(model, X);

% FGPLVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
%
%	Description:
%
%	[MU, SIGMA] = FGPLVMPOSTERIORMEANVAR(MODEL, X) returns the posterior
%	mean and variance for a given set of points.
%	 Returns:
%	  MU - the mean of the posterior distribution.
%	  SIGMA - the variances of the posterior distributions.
%	 Arguments:
%	  MODEL - the model for which the posterior will be computed.
%	  X - the input positions for which the posterior will be computed.
%	
%
%	See also
%	GPPOSTERIORMEANVAR, FGPLVMCREATE


%	Copyright (c) 2005, 2006 Neil D. Lawrence


[mu, varsigma] = gpPosteriorMeanVar(model, X);