function L = probitNoiseLikelihood(noise, mu, varsigma, y)


% PROBITNOISELIKELIHOOD Likelihood of the data under the PROBIT noise model.
% FORMAT
% DESC returns the likelihood of a data set under the  probit based classification noise model.
% ARG noise : the noise structure for which the likelihood is required.
% ARG mu : input mean locations for the likelihood.
% ARG varSigma : input variance locations for the likelihood.
% ARG y : target locations for the likelihood.
%
% SEEALSO : probitNoiseParamInit, probitNoiseLogLikelihood, noiseLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(y, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
L = cumGaussian((y.*mu)./(sqrt(noise.sigma2+varsigma)));
