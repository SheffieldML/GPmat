function [m, beta] = ngaussNoiseSites(noise, g, nu, mu, varSigma, y)

% NGAUSSNOISESITES Site updates for noiseless Gaussian noise model.

% SHEFFIELDML

% SHEFFIELDML


[m, beta] = gaussianNoiseSites(noise, g, nu, mu, varSigma, y);