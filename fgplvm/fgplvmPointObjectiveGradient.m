function f = fgplvmPointObjectiveGradient(x, model, y)

% FGPLVMPOINTOBJECTIVEGRADIENT Wrapper function for objective and gradient of a single point in latent space and the output location..
% FORMAT
% DESC provides a wrapper function for the negative log probability
% of a given data point under the posterior distribution of the
% Gaussian process induced by the training data. Also returns the
% gradient of the negative log probability with respect to the
% given latent point.
% ARG x : location in input space for the point.
% ARG model : the model structure for which the negative log
% probability of the given data under the posterior is to be computed.
% ARG y : the location in data space for the point.
% RETURN f : the negative of the log probability of the given data
% point under the posterior distribution induced by the training data.
% RETURN g : the gradient of the log probability with respect to
% the given latent point.
% 
% SEEALSO : fgplvmCreate, fgplvmPointLogLikelihood,
% fgplvmOptimisePoint, fgplvmObjective, fgplvmGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

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

