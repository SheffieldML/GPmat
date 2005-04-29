function [dlnZ_dmu, dlnZ_dvs] = noiseGradVals(noise, mu, varsigma, y)

% NOISEGRADVALS Gradient of noise model wrt mu and varsigma.

% NOISE

[dlnZ_dmu, dlnZ_dvs] = feval([noise.type 'NoiseGradVals'], noise, mu, ...
                             varsigma, y);
%/~
if any(isnan(dlnZ_dmu))
  warning('Gradient of mu is NaN')
end
if any(isnan(dlnZ_dvs))
  warning('Gradient of varsigma is NaN')
end
%~/