function g = negNoiseGradientParam(params, model, prior)

% NEGNOISEGRADIENTPARAM Wrapper function for calling noise gradients.

% NOISE

% NOISE


model.noise = noiseExpandParam(model.noise, params);
g = - feval([model.noise.type 'GradientParam'], model);

if prior
  g =g +params;
end