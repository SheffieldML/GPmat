function model = gaussianNoiseParamInit(model)

% GAUSSIANNOISEPARAMINIT Gaussian noise model's parameter initialisation.

% IVM

model.noise.sigma2 = 1;
model.noise.bias = mean(model.y);