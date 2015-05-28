function L = mgaussianNoiseLogLikelihood(noise, mu, varsigma, y)


% MGAUSSIANNOISELOGLIKELIHOOD Log likelihood of the data under the MGAUSSIAN noise model.
% FORMAT
% DESC returns the log likelihood of a data set under the  multiple output Gaussian noise model.
% ARG noise : the noise structure for which the log likelihood is required.
% ARG mu : input mean locations for the log likelihood.
% ARG varSigma : input variance locations for the log likelihood.
% ARG y : target locations for the log likelihood.
%
% SEEALSO : mgaussianNoiseParamInit, mgaussianNoiseLikelihood, noiseLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% NOISE


N = size(mu, 1);
D = size(mu, 2);
for i = 1:D
  varsigma(:, i) = varsigma(:, i) + noise.sigma2(i);
  mu(:, i) = mu(:, i) + noise.bias(i);
end
arg = (y - mu);
arg = arg.*arg./varsigma;

% Remove unlabelled data from likelihood.
arg = arg(:);
unlabelled = find(isnan(arg));
arg(unlabelled) = [];
varsigma(unlabelled) = [];
L = - 0.5*sum(sum(log(varsigma))) ...
    - 0.5*sum(sum(arg)) ...
    - 0.5*N*D*log(2*pi);
