function f = fgplvmPointObjective(x, model, y)

% FGPLVMPOINTOBJECTIVE Wrapper function for objective of a single point in latent space and the output location..
% FORMAT
% DESC provides a wrapper function for the negative log probability
% of a given data point under the posterior distribution of the
% Gaussian process induced by the training data.
% ARG x : location in input space for the point.
% ARG model : the model structure for which the negative log
% probability of the given data under the posterior is to be computed.
% ARG y : the location in data space for the point.
% RETURN f : the negative of the log probability of the given data
% point under the posterior distribution induced by the training data.
% 
% SEEALSO : fgplvmCreate, fgplvmPointLogLikelihood, fgplvmOptimisePoint
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

f = - fgplvmPointLogLikelihood(model, x, y);
