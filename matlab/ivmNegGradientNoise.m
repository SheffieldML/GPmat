function g = ivmNegGradientNoise(params, model)

% IVMNEGGRADIENTNOISE Wrapper function for calling noise param gradients.

% IVM

model.noise = noiseExpandParam(model.noise, params);
g = - noiseGradientParam(model.noise, model.mu, model.varSigma, model.y);
