function params = heavisideNoiseExtractParam(noise)

% HEAVISIDENOISEEXTRACTPARAM Extract parameters from heaviside noise model.

% IVM

params = [invSigmoid(noise.eta) noise.bias];

