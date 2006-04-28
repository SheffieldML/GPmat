function y = gpOut(model, x);

% GPOUT Evaluate the output of an Gaussian process model.
%
% y = gpOut(model, x);
%

% Copyright (c) 2006 Neil D. Lawrence
% gpOut.m version 1.1



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
