function noise = noiseParamInit(noise, y)

% NOISEPARAMINIT Noise model's parameter initialisation.

% IVM

noise = feval([noise.type 'NoiseParamInit'], noise, y);