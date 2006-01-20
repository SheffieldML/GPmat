function [mu, varsigma] = fgplvmDynamicsPosteriorMeanVar(model, X);

% FGPLVMDYNAMICSPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.

% FGPLVM

if ~isfield(model, 'alpha')
  if model.dynamics.diff
    Y = model.X(2:end, :) - model.X(1:end-1, :);
  else
    Y = model.X(2:end, :);
  end
  model.dynamics = gpComputeAlpha(model.dynamics, Y);
end
[mu, varsigma] = gpPosteriorMeanVar(model.dynamics, X);
