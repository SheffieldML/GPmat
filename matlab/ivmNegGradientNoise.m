function g = ivmNegGradientNoise(params, model)

% IVMNEGGRADIENTNOISE Wrapper function for calling noise param gradients.
% FORMAT
% DESC is a wrapper function that returns the negative gradients of
% the log likelihood with respect to the noise parameters for
% optimisation in the NETLAB style.
% ARG params : the parameters where the gradients are to be
% computed.
% ARG model : the model for which the gradients are to be computed.
% RETURN g : the negative gradients of the log likelihood with respect to
% the noise parameters.
%
% SEEALSO : noiseGradientParam, noiseExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005
% 

% IVM

model.noise = noiseExpandParam(model.noise, params);
g = - noiseGradientParam(model.noise, model.mu, model.varSigma, model.y);
