function f = fgplvmSequenceObjective(xvec, model, Y, varargin)

% FGPLVMSEQUENCEOBJECTIVE Wrapper function for objective of a single sequence in latent space and the corresponding output sequence.
%
%	Description:
%
%	F = FGPLVMSEQUENCEOBJECTIVE(VECX, MODEL, Y, ...) provides a wrapper
%	function for the negative log probability of a given data sequence
%	under the posterior distribution of the Gaussian process induced by
%	the training data..
%	 Returns:
%	  F - the negative of the log probability of the given data sequence
%	   under the posterior distribution induced by the training data.
%	 Arguments:
%	  VECX - time ordered locations in input space for the sequence
%	   placed as a vector using the matlab X(:)' notation.
%	  MODEL - the model structure for which the negative log probability
%	   of the given data under the posterior is to be computed.
%	  Y - time ordered locations in data spaces for the sequence.
%	  ... - optional additional arguments to be passed to the model
%	   sequence log likelihood.
%	
%
%	See also
%	FGPLVMCREATE, FGPLVMSEQUENCELOGLIKELIHOOD, FGPLVMOPTIMISESEQUENCE


%	Copyright (c) 2006 Neil D. Lawrence


X = reshape(xvec, size(Y, 1), model.q);
f = - fgplvmSequenceLogLikelihood(model, X, Y, varargin{:});
