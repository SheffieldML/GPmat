function e = negIvmLogLikelihood(params, model, prior)

% NEGIVMLOGLIKELIHOOD Wrapper function for calling ivm likelihoods.


model.noise = noiseExpandParam(model.noise, params);
e = - ivmLogLikelihood(model);

if prior
  e =e +0.5*params*params';
end