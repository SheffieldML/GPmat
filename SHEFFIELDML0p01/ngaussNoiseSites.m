function [m, beta] = ngaussNoiseSites(noise, g, nu, mu, varSigma, y)

% NGAUSSNOISESITES Site updates for noiseless Gaussian noise model.
%
%	Description:
%	[m, beta] = ngaussNoiseSites(noise, g, nu, mu, varSigma, y)
%



[m, beta] = gaussianNoiseSites(noise, g, nu, mu, varSigma, y);