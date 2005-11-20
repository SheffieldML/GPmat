function y = ivmOut(model, x);

% IVMOUT Evaluate the output of an ivm model.

% IVM 

if nargin < 2
  % This implies evaluate for the training data.
  mu = model.mu;
  varsigma = model.varSigma;
else
  [mu, varsigma] = ivmPosteriorMeanVar(model, x);
end

y = noiseOut(model.noise, mu, varsigma);
