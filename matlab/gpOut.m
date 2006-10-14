function y = gpOut(model, x);

% GPOUT Evaluate the output of an Gaussian process model.

% FGPLVM 

if nargin < 2
  % This implies evaluate for the training data.
  mu = model.mu;
  varsigma = model.varSigma;
else
  [mu, varsigma] = gpPosteriorMeanVar(model, x);
end

if isfield(model, 'noise')
  y = noiseOut(model.noise, mu, varsigma);
end
