function g = fgplvmSequenceGradient(xvec, model, Y, varargin)

% FGPLVMSEQUENCEGRADIENT Wrapper function for gradient of a latent sequence.
%
%	Description:
%
%	G = FGPLVMSEQUENCEGRADIENT(VECX, MODEL, Y, ...) is a wrapper
%	function for the gradient of the log probability of a sequence under
%	the posterior distribution induced by the training data with respect
%	to the latent positions.  The GP-LVM model is one that is already
%	trained with a specific data set.
%	 Returns:
%	  G - the gradient of the log likelihood with respect to the latent
%	   positions.
%	 Arguments:
%	  VECX - the positions in the latent space that are being optimised
%	   as a row vector using the matlab X(:)' notation.
%	  MODEL - the trained GP-LVM model that is being optimised.
%	  Y - the position in data space for which the latent sequence is
%	   being optimised.
%	  ... - optional additional arguments to be passed to the model
%	   sequence log likelihood gradient.
%	
%
%	See also
%	FGPLVMSEQUENCELOGLIKEGRADIENT, FGPLVMOPTIMISESEQUENCE


%	Copyright (c) 2006 % COPYRIGHT Neil D. Lawrence


X = reshape(xvec, size(Y, 1), model.q);
g = - fgplvmSequenceLogLikeGradient(model, X, Y, varargin{:});

g = g(:)';