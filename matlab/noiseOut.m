function y = noiseOut(noise, mu, varsigma);

% NOISEOUT Give the output of the noise model given the mean and variance.

% NOISE

fhandle = str2func([noise.type 'NoiseOut']);
y = fhandle(noise, mu, varsigma);
