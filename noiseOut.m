function y = noiseOut(noise, mu, varsigma);

% NOISEOUT Give the output of the noise model given the mean and variance.
% FORMAT
% DESC computes the ouptut for the given
% noise given input mean and variances.
% ARG noise : the noise structure for which the output is computed.
% ARG mu : the input mean values.
% ARG varSigma : the input variance values.
% RETURN y : the output from the noise model.
%
% SEEALSO : noiseParamInit, noiseCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

fhandle = str2func([noise.type 'NoiseOut']);
y = fhandle(noise, mu, varsigma);
