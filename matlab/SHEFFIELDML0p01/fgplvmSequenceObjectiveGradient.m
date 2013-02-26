function [f, g] = fgplvmSequenceObjectiveGradient(xvec, model, Y)

% FGPLVMSEQUENCEOBJECTIVEGRADIENT Wrapper function for objective
%
%	Description:
%	and gradient of a single sequence in latent space and the corresponding output sequence.
%
%	[F, G] = FGPLVMSEQUENCEOBJECTIVEGRADIENT(VECX, MODEL, Y) provides a
%	wrapper function for the negative log probability of a given data
%	sequence under the posterior distribution of the Gaussian process
%	induced by the training data..
%	 Returns:
%	  F - the negative of the log probability of the given data sequence
%	   under the posterior distribution induced by the training data.
%	  G - the gradient of the negative of the returned log probability
%	   with respect to the latent sequence.
%	 Arguments:
%	  VECX - time ordered locations in input space for the sequence
%	   placed as a vector using the matlab X(:)' notation.
%	  MODEL - the model structure for which the negative log probability
%	   of the given data under the posterior is to be computed.
%	  Y - time ordered locations in data spaces for the sequence.
%	
%
%	See also
%	FGPLVMCREATE, FGPLVMSEQUENCELOGLIKELIHOOD, FGPLVMOPTIMISESEQUENCE


%	Copyright (c) 2006 Neil D. Lawrence


% Check how the optimiser has given the parameters
if size(xvec, 1) > size(xvec, 2)
  % As a column vector ... transpose everything.
  transpose = true;
  X = reshape(xvec', size(Y, 1), model.q);
else
  transpose = false;
  X = reshape(xvec, size(Y, 1), model.q);
end
f = - fgplvmSequenceLogLikelihood(model, X, Y);

if nargout > 1
  g = - fgplvmSequenceLogLikeGradient(model, X, Y);
  g = g(:)';  
end
if transpose
  g = g';
end