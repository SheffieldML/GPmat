function [varsigma] = ppcaPosteriorVar(model, X);

% PPCAPOSTERIORVAR Mean and variances of the posterior at points given by X.
%
%	Description:
%
%	SIGMA = PPCAPOSTERIORVAR(MODEL, X) returns the posterior mean and
%	variance for a given set of points.
%	 Returns:
%	  SIGMA - the variances of the posterior distributions.
%	 Arguments:
%	  MODEL - the model for which the posterior will be computed.
%	  X - the input positions for which the posterior will be computed.
%	
%
%	See also
%	PPCACREATE, PPCAPOSTERIORMEANVAR


%	Copyright (c) 2008 Neil D. Lawrence


varsigma = repmat(1/model.beta, size(X, 1), 1);
