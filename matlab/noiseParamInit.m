function model = noiseParamInit(model)

% NOISEPARAMINIT Noise model's parameter initialisation.

% IVM

model = feval([model.noise.type 'NoiseParamInit'], model);