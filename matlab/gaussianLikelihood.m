function L = gaussianLikelihood(noise, mu, varsigma, y)

% GAUSSIANLIKELIHOOD Likelihood of data under Gaussian noise model.

% NOISE

N = size(y, 1);
D = size(y, 2);
varsigma = varsigma + noise.sigma2;
for i = 1:D
  mu(:, i) = mu(:, i) + noise.bias(i);
end
arg = (mu - y)./sqrt(varsigma);
L = (2*pi*varsigma).^(-1/2).*exp( - .5*arg.*arg);
