function g = fgplvmGradient(params, model)

% FGPLVMGRADIENT GP-LVM gradient wrapper.
% FORMAT
% DESC is a wrapper function for the gradient of the negative log
% likelihood of an GP-LVM model with respect to the latent postions
% and parameters.
% ARG params : vector of parameters and latent postions where the
% gradient is to be evaluated.
% ARG model : the model structure into which the latent positions
% and the parameters will be placed.
% RETURN g : the gradient of the negative log likelihood with
% respect to the latent positions and the parameters at the given
% point.
% 
% SEEALSO : fgplvmLogLikeGradients, fgplvmExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2005

% FGPLVM

model = fgplvmExpandParam(model, params);
g = - fgplvmLogLikeGradients(model);
