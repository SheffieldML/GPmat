function L = ivmLikelihoods(model, x, y);

% IVMLOGLIKELIHOODS Return the likelihood for each point for the IVM.

if nargin < 3
  % This implies evaluate for the traing data.
  mu = model.mu;
  varsigma = model.varSigma;
  y = model.y;
else
  [mu, varsigma] = ivmPosteriorMeanVar(model, x);
end

L = noiseLikelihood(model.noise, mu, varsigma, y);
