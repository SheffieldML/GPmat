function L = noiseLogLikelihood(noise, mu, varsigma, y);

% NOISELOGLIKELIHOOD Return the log-likelihood under the noise model.

L = feval([noise.type 'LogLikelihood'], noise, mu, varsigma, y);
