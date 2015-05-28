function L = noiseExpectationLogLikelihood(noise, mu, varsigma, y);

% NOISEEXPECTATIONLOGLIKELIHOOD Return the expectation of the log likelihood.
% FORMAT
% DESC returns the expectation of the log likelihood for a gven noise model.
% ARG noise : the noise structure for which the expectation of the log
% likelihood is required.
% ARG mu : input mean locations for the likelihood.
% ARG varSigma : input variance locations for the likelihood.
% ARG y : target locations for the likelihood.
%
% SEEALSO : noiseParamInit, noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2007

% NOISE

fhandle = str2func([noise.type 'NoiseExpectationLogLikelihood']);
L = fhandle(noise, mu, varsigma, y);


