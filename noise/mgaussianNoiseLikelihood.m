function L = mgaussianNoiseLikelihood(noise, mu, varsigma, y)


% MGAUSSIANNOISELIKELIHOOD Likelihood of the data under the MGAUSSIAN noise model.
% FORMAT
% DESC returns the likelihood of a data set under the  multiple output Gaussian noise model.
% ARG noise : the noise structure for which the likelihood is required.
% ARG mu : input mean locations for the likelihood.
% ARG varSigma : input variance locations for the likelihood.
% ARG y : target locations for the likelihood.
%
% SEEALSO : mgaussianNoiseParamInit, mgaussianNoiseLogLikelihood, noiseLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


N = size(y, 1);
D = size(y, 2);
for i = 1:D
  varsigma(:, i) = varsigma(:, i) + noise.sigma2(i);
  mu(:, i) = mu(:, i) + noise.bias(i);
end
arg = (mu - y)./sqrt(varsigma);

L = (2*pi*varsigma).^(-1/2).*exp( - .5*arg.*arg);

% Set likelihood of unlabelled points to 1.
L(find(isnan(y)) = 1;