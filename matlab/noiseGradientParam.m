function g = noiseGradientParam(noise, mu, varsigma, y)

% NOISEGRADIENTPARAM Gradient of the noise model's parameters.

% IVM

g = feval([noise.type 'NoiseGradientParam'], noise, mu, varsigma, y);