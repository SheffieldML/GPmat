function g = fgplvmSequenceGradient(xvec, model, Y, varargin)

% FGPLVMSEQUENCEGRADIENT Wrapper function for gradient of a latent sequence.
% FORMAT
% DESC is a wrapper function for the gradient of the log probability
% of a sequence under the posterior distribution induced by the
% training data with respect to the latent positions.  The GP-LVM
% model is one that is already trained with a specific data set.
% ARG vecx : the positions in the latent space that are being
% optimised as a row vector using the matlab X(:)' notation.
% ARG model : the trained GP-LVM model that is being optimised.
% ARG Y : the position in data space for which the latent sequence is
% being optimised.
% ARG P1, P2, P3 ... : optional additional arguments to be passed
% to the model sequence log likelihood gradient.
% RETURN g : the gradient of the log likelihood with respect to the
% latent positions.
%
% SEEALSO : fgplvmSequenceLogLikeGradient, fgplvmOptimiseSequence
%
% COPYRIGHT Neil D. Lawrence, 2006

% FGPLVM

X = reshape(xvec, size(Y, 1), model.q);
g = - fgplvmSequenceLogLikeGradient(model, X, Y, varargin{:});

g = g(:)';
