function e = negNoiseLogLikelihood(params, model, prior)

% NEGNOISELOGLIKELIHOOD Wrapper function for calling noise likelihoods.

% NOISE

model.noise = noiseExpandParam(model.noise, params);
fhandle = str2func([model.noise.type 'NoiseLogLikelihood']);
e = - fhandle([], [], model);

if prior
  e =e +0.5*params*params';
end