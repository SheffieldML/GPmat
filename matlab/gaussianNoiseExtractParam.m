function params = gaussianNoiseExtractParam(noise)

% GAUSSIANNOISEEXTRACTPARAM Extract parameters from Gaussian noise model.

% IVM

params = [noise.bias log(noise.sigma2)];