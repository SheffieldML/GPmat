function L = noiseLogLikelihood(noise, mu, varsigma, y);

% NOISELOGLIKELIHOOD Return the log-likelihood under the noise model.

% IVM

L = feval([noise.type 'LogLikelihood'], noise, mu, varsigma, y);
