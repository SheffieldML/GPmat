function L = mgaussianLogLikelihood(noise, mu, varsigma, y)

% MGAUSSIANLOGLIKELIHOOD Log-likelihood of data under Variable variance Gaussian noise model.

% NOISE

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
