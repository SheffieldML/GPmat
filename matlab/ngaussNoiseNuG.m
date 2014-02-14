function [g, nu] = ngaussNoiseNuG(noise, mu, varSigma, y)

% NGAUSSNOISENUG Update nu and g parameters associated with noiseless Gaussian noise model.

% GPMAT

% GPMAT


[g, nu] = gaussianNoiseNuG(noise, mu, varSigma, y);
