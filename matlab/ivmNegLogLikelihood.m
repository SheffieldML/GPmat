function e = ivmNegLogLikelihood(params, model)

% IVMNEGLOGLIKELIHOOD Wrapper function for calling ivm likelihoods.

% IVM

model.noise = noiseExpandParam(model.noise, params);
e = - ivmLogLikelihood(model);
