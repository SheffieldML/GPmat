function varsigma = fgplvmPosteriorVar(model, X);

% FGPLVMPOSTERIORVAR Variances of the posterior at points given by X.
%
%	Description:
%
%	SIGMA = FGPLVMPOSTERIORVAR(MODEL, X) returns the posterior mean and
%	variance for a given set of points.
%	 Returns:
%	  SIGMA - the variances of the posterior distributions.
%	 Arguments:
%	  MODEL - the model for which the posterior will be computed.
%	  X - the input positions for which the posterior will be computed.
%	
%
%	See also
%	GPPOSTERIORVAR, FGPLVMCREATE


%	Copyright (c) 2009 Neil D. Lawrence


varsigma = gpPosteriorVar(model, X);