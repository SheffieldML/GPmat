function [g, nu] = gaussianNoiseNuG(noise, mu, varSigma, y)

% GAUSSIANNOISENUG Compute nu and g for GAUSSIAN noise model.
% FORMAT
% DESC computes the values nu and g for the Gaussian noise given the mean and variance inputs as well as the output of the noise model.
% ARG noise : the noise structure for which the nu and g are computed.
% ARG mu : input mean to the noise model.
% ARG varSigma : input variance to the noise model.
% ARG y : target output for the noise model.
% RETURN g : the vector g, which is the gradient of log Z with respect to the input mean.
% ARG y : target output for the noise model.
% RETURN nu : the vector nu, see equation 10 of "Extensions of the Informative Vector Machine".
%
% SEEALSO : gaussianNoiseParamInit, noiseUpdateNuG, noiseCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(y, 2);
nu = 1./(noise.sigma2+varSigma);
g = zeros(size(nu));
for i = 1:D
  g(:, i) = y(:, i) - mu(:, i) - ...
      noise.bias(i);
end
g = g.*nu;
