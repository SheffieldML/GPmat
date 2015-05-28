function L = gaussianLogLikelihood(noise, mu, varsigma, y)

% GAUSSIANLOGLIKELIHOOD Log-likelihood of data under Gaussian noise model.

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

