function f = fgplvmPointObjectiveGradient(x, model, y)

% FGPLVMPOINTOBJECTIVEGRADIENT Wrapper function for objective and gradient of a single point in latent space and the output location..
%
%	Description:
%
%	[F, G] = FGPLVMPOINTOBJECTIVEGRADIENT(X, MODEL, Y) provides a
%	wrapper function for the negative log probability of a given data
%	point under the posterior distribution of the Gaussian process
%	induced by the training data. Also returns the gradient of the
%	negative log probability with respect to the given latent point.
%	 Returns:
%	  F - the negative of the log probability of the given data point
%	   under the posterior distribution induced by the training data.
%	  G - the gradient of the log probability with respect to the given
%	   latent point.
%	 Arguments:
%	  X - location in input space for the point.
%	  MODEL - the model structure for which the negative log probability
%	   of the given data under the posterior is to be computed.
%	  Y - the location in data space for the point.
%	fgplvmOptimisePoint, fgplvmObjective, fgplvmGradient
%	
%
%	See also
%	FGPLVMCREATE, FGPLVMPOINTLOGLIKELIHOOD, 


%	Copyright (c) 2006 Neil D. Lawrence


% Check how the optimiser has given the parameters
if size(xvec, 1) > size(xvec, 2)
  % As a column vector ... transpose everything.
  transpose = true;
  x = x';
else
  transpose = false;
end
f = - fgplvmPointLogLikelihood(model, x, y);

if nargout > 1
  g = - fgplvmPointLogLikeGradient(model, x, y);
end
if transpose
  g = g';
end

