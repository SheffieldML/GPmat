function g = negNoiseGradientParam(params, model, prior)

% NEGNOISEGRADIENTPARAM Wrapper function for calling noise gradients.

% IVM

model.noise = noiseExpandParam(params, model.noise);
g = - feval([model.noise.type 'GradientParam'], model);

if prior
  g =g +params;
end