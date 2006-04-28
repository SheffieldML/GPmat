function y = gpDynamicsSamp(model, X);

% GPDYNAMICSSAMP Sample from the dynamics for a given input.
%
% y = gpDynamicsSamp(model, X);
%

% Copyright (c) 2006 Neil D. Lawrence
% gpDynamicsSamp.m version 1.1



[mu, var] = gpPosteriorMeanVar(model, X);
y = gsamp(mu, diag(var), 1);
if isfield(model, 'diff') & model.diff
  y = X + y;
end