function [dlnZ_dmu, dlnZ_dvs] = noiseGradVals(noise, mu, varsigma, y)

% NOISEGRADVALS Gradient of noise model wrt mu and varsigma.

% IVM

[dlnZ_dmu, dlnZ_dvs] = feval([noise.type 'NoiseGradVals'], noise, mu, varsigma, y);