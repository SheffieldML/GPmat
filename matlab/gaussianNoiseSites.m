function [m, beta] = gaussianNoiseSites(noise, g, nu, mu, varSigma, y)

% GAUSSIANNOISESITES Update the site parameters for the GAUSSIAN noise mode.
% FORMAT
% DESC updates the site parameters for the Gaussian
% noise model. 
% ARG noise : the noise structure for which the site parameters are to
% be updated. 
% ARG g : values of g as retuned by gaussianNoiseNuG.
% ARG nu : values of nu as retuned by gaussianNoiseNuG.
% ARG mu : the mean value of the Gaussian input to the noise structure.
% ARG varSigma : the variance of the Gaussian input to the noise structure.
% ARG y : the target value.
% RETURN m : the site mean parameters.
% RETURN beta : the site precision parameters.
%
% SEEALSO : gaussianNoiseParamInit, noiseUpdateSites
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

N = size(y, 1);
D = length(noise.bias);
beta = zeros(N, D);
for i = 1:size(y, 2)
  m(:, i) = y(:, i) - noise.bias(i);
end
beta = repmat(1./noise.sigma2, N, D);