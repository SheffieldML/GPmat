function L = noiseLikelihood(noise, mu, varsigma, y);

% NOISELIKELIHOOD Return the likelihood for each point under the noise model.

% NOISE

fhandle = str2func([noise.type 'Likelihood']);
L = fhandle(noise, mu, varsigma, y);
