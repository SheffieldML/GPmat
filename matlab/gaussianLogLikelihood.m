function L = gaussianLogLikelihood(X, Y, model, prior)

% GAUSSIANLOGLIKELIHOOD Log-likelihood of data under Gaussian noise model.

% IVM

if isempty(X)
  mu = model.mu;
  varSigma = model.varSigma;
  Y = model.y;
else
  [mu, varSigma] = ivmPosteriorMeanVar(X, model);
end
N = size(Y, 1);
D = size(Y, 2);
varSigma = varSigma + model.noise.sigma2;
for i = 1:D
  mu(:, i) = mu(:, i) + model.noise.bias(i);
end
arg = (Y - mu);
arg = arg.*arg./varSigma;

L = - 0.5*sum(sum(log(varSigma))) ...
    - 0.5*sum(sum(arg)) ...
    - 0.5*N*D*log(2*pi);

