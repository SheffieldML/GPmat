function [mu, varsigma] = fgplvmDynamicsPosteriorMeanVar(model, X);

% FGPLVMDYNAMICSPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
%
% [mu, varsigma] = fgplvmDynamicsPosteriorMeanVar(model, X);
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmDynamicsPosteriorMeanVar.m version 1.2



if ~isfield(model, 'alpha')
  if model.dynamics.diff
    Y = model.X(2:end, :) - model.X(1:end-1, :);
  else
    Y = model.X(2:end, :);
  end
  model.dynamics = gpComputeAlpha(model.dynamics, Y);
end
[mu, varsigma] = gpPosteriorMeanVar(model.dynamics, X);
