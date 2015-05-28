function L = noiseLikelihood(noise, mu, varsigma, y);

% NOISELIKELIHOOD Return the likelihood for each point under the noise model.
% FORMAT
% DESC returns the likelihoods for data points under the given noise model.
% ARG noise : the noise structure for which the likelihood is required.
% ARG mu : input mean locations for the likelihood.
% ARG varSigma : input variance locations for the likelihood.
% ARG y : target locations for the likelihood.
%
% SEEALSO : noiseParamInit, noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE

fhandle = str2func([noise.type 'NoiseLikelihood']);
L = fhandle(noise, mu, varsigma, y);
