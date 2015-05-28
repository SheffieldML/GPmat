function L = gaussianNoiseExpectationLogLikelihood(noise, mu, varsigma, y)

% GAUSSIANNOISEEXPECTATIONLOGLIKELIHOOD Likelihood of the data under the GAUSSIAN noise model.
% FORMAT
% DESC returns the likelihoods for data points under the  Gaussian noise model.
% ARG noise : the noise structure for which the likelihood is required.
% ARG mu : input mean locations for the likelihood.
% ARG varSigma : input variance locations for the likelihood.
% ARG y : target locations for the likelihood.
%
% SEEALSO : gaussianNoiseParamInit, gaussianNoiseLogLikelihood, noiseLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005



N = size(y, 1);
D = size(y, 2);
varsigma = varsigma + noise.sigma2;
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
arg = (mu - y)./sqrt(varsigma);
L = (2*pi*varsigma).^(-1/2).*exp( - .5*arg.*arg);
