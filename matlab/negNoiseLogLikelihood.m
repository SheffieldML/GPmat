function e = negNoiseLogLikelihood(params, model, prior)

% NEGNOISELOGLIKELIHOOD Wrapper function for calling noise likelihoods.

% NOISE

% NOISE


model.noise = noiseExpandParam(model.noise, params);
e = - feval([model.noise.type 'LogLikelihood'], [], [], model);

if prior
  e =e +0.5*params*params';
end