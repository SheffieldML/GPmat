function L = noiseLogLikelihood(x, y, model);

% NOISELOGLIKELIHOOD Return the log-likelihood under the noise model.

L = feval([model.noise.type 'LogLikelihood'], x, y, model);