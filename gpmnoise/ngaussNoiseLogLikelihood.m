function L = ngaussNoiseLogLikelihood(noise, mu, varsigma, y)


% NGAUSSNOISELOGLIKELIHOOD Log likelihood of the data under the NGAUSS noise model.
% FORMAT
% DESC returns the log likelihood of a data set under the  noiseless Gaussian noise model.
% ARG noise : the noise structure for which the log likelihood is required.
% ARG mu : input mean locations for the log likelihood.
% ARG varSigma : input variance locations for the log likelihood.
% ARG y : target locations for the log likelihood.
%
% SEEALSO : ngaussNoiseParamInit, ngaussNoiseLikelihood, noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


L = gaussianNoiseLogLikelihood(noise, mu, varsigma, y);
