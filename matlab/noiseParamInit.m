function noise = noiseParamInit(noise, y)

% NOISEPARAMINIT Noise model's parameter initialisation.

% IVM

if nargin > 1
  noise = feval([noise.type 'NoiseParamInit'], noise, y);
else
  noise = feval([noise.type 'NoiseParamInit'], noise);
end
