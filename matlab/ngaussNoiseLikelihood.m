function L = ngaussNoiseLikelihood(noise, mu, varsigma, y)


% NGAUSSNOISELIKELIHOOD Likelihood of the data under the NGAUSS noise model.
% FORMAT
% DESC returns the likelihood of a data set under the  noiseless Gaussian noise model.
% ARG noise : the noise structure for which the likelihood is required.
% ARG mu : input mean locations for the likelihood.
% ARG varSigma : input variance locations for the likelihood.
% ARG y : target locations for the likelihood.
%
% SEEALSO : ngaussNoiseParamInit, ngaussNoiseLogLikelihood, noiseLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


L = gaussianNoiseLikelihood(noise, mu, varsigma, y);