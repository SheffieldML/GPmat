function [g, nu] = ngaussNoiseNuG(noise, mu, varSigma, y)

% NGAUSSNOISENUG Update nu and g parameters associated with noiseless Gaussian noise model.

% SHEFFIELDML

% SHEFFIELDML


[g, nu] = gaussianNoiseNuG(noise, mu, varSigma, y);
