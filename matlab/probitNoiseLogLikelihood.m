function L = probitNoiseLogLikelihood(noise, mu, varsigma, y)


% PROBITNOISELOGLIKELIHOOD Log likelihood of the data under the PROBIT noise model.
% FORMAT
% DESC returns the log likelihood of a data set under the  probit based classification noise model.
% ARG noise : the noise structure for which the log likelihood is required.
% ARG mu : input mean locations for the log likelihood.
% ARG varSigma : input variance locations for the log likelihood.
% ARG y : target locations for the log likelihood.
%
% SEEALSO : probitNoiseParamInit, probitNoiseLikelihood, noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


D = size(y, 2);
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end

L = sum(sum(lnCumGaussian((y.*mu)./(sqrt(noise.sigma2+varsigma)))));
