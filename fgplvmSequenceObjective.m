function f = fgplvmSequenceObjective(xvec, model, Y, varargin)

% FGPLVMSEQUENCEOBJECTIVE Wrapper function for objective of a single sequence in latent space and the corresponding output sequence.
% FORMAT
% DESC provides a wrapper function for the negative log probability
% of a given data sequence under the posterior distribution of the
% Gaussian process induced by the training data..
% ARG vecx : time ordered locations in input space for the sequence
% placed as a vector using the matlab X(:)' notation.
% ARG model : the model structure for which the negative log
% probability of the given data under the posterior is to be computed.
% ARG Y : time ordered locations in data spaces for the sequence.
% ARG P1, P2, P3 ... : optional additional arguments to be passed
% to the model sequence log likelihood.
% RETURN f : the negative of the log probability of the given data
% sequence under the posterior distribution induced by the training data.
% 
% SEEALSO : fgplvmCreate, fgplvmSequenceLogLikelihood, fgplvmOptimiseSequence
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

X = reshape(xvec, size(Y, 1), model.q);
f = - fgplvmSequenceLogLikelihood(model, X, Y, varargin{:});
