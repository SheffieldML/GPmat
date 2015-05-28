function L = mgaussianLikelihood(noise, mu, varsigma, y)

% MGAUSSIANLIKELIHOOD Likelihood of data under Variable variance Gaussian noise model.

% NOISE

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