function g = fgplvmPointGradient(x, model, y)

% FGPLVMPOINTGRADIENT Wrapper function for gradient of a single point.
% FORMAT
% DESC is a wrapper function for the gradient of the log likelihood
% with respect to a point in the latent space. The GP-LVM
% model is one that is assumed to have already been trained.
% ARG x : the position in the latent space that is being optimised.
% ARG model : the trained GP-LVM model that is being optimised.
% ARG y : the position in data space for which the latent point is
% being optimised.
% RETURN g : the gradient of the log likelihood with respect to the
% latent position.
%
% SEEALSO : fgplvmPointLogLikeGradient, fgplvmOptimisePoint
%
% COPYRIGHT Neil D. Lawrence, 2005, 2006

% FGPLVM

g = - fgplvmPointLogLikeGradient(model, x, y);
