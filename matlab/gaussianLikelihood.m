function L = gaussianLikelihood(X, Y, model)

% GAUSSIANLIKELIHOOD Likelihood of data under Gaussian noise model.

% IVM

N = size(Y, 1);
D = size(Y, 2);
[mu, varsigma] = ivmPosteriorMeanVar(X, model);
varsigma = varsigma + model.noise.sigma2;
for i = 1:D
  mu(:, i) = mu(:, i) + model.noise.bias(i);
end
arg = (mu - Y)./sqrt(varsigma);
L = (2*pi*varsigma).^(-1/2).*exp( - .5*arg.*arg);
