function noise = gaussianNoiseParamInit(noise, y)

% GAUSSIANNOISEPARAMINIT Gaussian noise model's parameter initialisation.

% IVM

noise.sigma2 = 1;
noise.bias = mean(y);
noise.nParams = 1 + size(y, 2);