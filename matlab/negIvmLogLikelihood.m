function e = negIvmLogLikelihood(params, model, prior)

% NEGIVMLOGLIKELIHOOD Wrapper function for calling ivm likelihoods.

% IVM

model.noise = noiseExpandParam(params, model.noise);
e = - ivmLogLikelihood(model);

if prior
  e =e +0.5*params*params';
end