function [params, names] = ngaussNoiseExtractParam(noise)

% NGAUSSNOISEEXTRACTPARAM Extract parameters from noiseless Gaussian noise model.

% IVM

params = [noise.bias];