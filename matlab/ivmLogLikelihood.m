function L = ivmLogLikelihood(model, x, y);

% IVMLOGLIKELIHOOD Return the log-likelihood for the IVM.

% IVM

if nargin < 3
  % This implies evaluate for the traing data.
  mu = model.mu;
  varsigma = model.varSigma;
  y = model.y;
else
  [mu, varsigma] = ivmPosteriorMeanVar(model, x);
end

L = noiseLogLikelihood(model.noise, mu, varsigma, y);
