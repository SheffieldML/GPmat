function f = fgplvmPointObjective(x, model, y)

% FGPLVMPOINTOBJECTIVE Wrapper function for objective of a single point in latent space and the output location..
%
%	Description:
%
%	F = FGPLVMPOINTOBJECTIVE(X, MODEL, Y) provides a wrapper function
%	for the negative log probability of a given data point under the
%	posterior distribution of the Gaussian process induced by the
%	training data.
%	 Returns:
%	  F - the negative of the log probability of the given data point
%	   under the posterior distribution induced by the training data.
%	 Arguments:
%	  X - location in input space for the point.
%	  MODEL - the model structure for which the negative log probability
%	   of the given data under the posterior is to be computed.
%	  Y - the location in data space for the point.
%	
%
%	See also
%	FGPLVMCREATE, FGPLVMPOINTLOGLIKELIHOOD, FGPLVMOPTIMISEPOINT


%	Copyright (c) 2006 Neil D. Lawrence


f = - fgplvmPointLogLikelihood(model, x, y);
