function g = negIvmGradientNoise(params, model, prior)

% NEGIVMGRADIENTNOISE Wrapper function for calling noise param gradients.

% IVM

model.noise = noiseExpandParam(model.noise, params);
g = - noiseGradientParam(model.noise, model.mu, model.varSigma, model.y);

if prior
  g =g +params;
end