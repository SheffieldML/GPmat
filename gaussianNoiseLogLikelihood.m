function L = gaussianNoiseLogLikelihood(noise, mu, varsigma, y)

% GAUSSIANNOISELOGLIKELIHOOD Log likelihood of the data under the GAUSSIAN noise model.
% FORMAT
% DESC returns the log likelihood of a data set under the  Gaussian noise model.
% ARG noise : the noise structure for which the log likelihood is required.
% ARG mu : input mean locations for the log likelihood.
% ARG varSigma : input variance locations for the log likelihood.
% ARG y : target locations for the log likelihood.
%
% SEEALSO : gaussianNoiseParamInit, gaussianNoiseLikelihood, noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


N = size(mu, 1);
D = size(mu, 2);
varsigma = varsigma + noise.sigma2;
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
arg = (y - mu);
arg = arg.*arg./varsigma;

L = - 0.5*sum(sum(log(varsigma))) ...
    - 0.5*sum(sum(arg)) ...
    - 0.5*N*D*log(2*pi);

