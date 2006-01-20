function y = gpDynamicsSamp(model, X);

% GPDYNAMICSSAMP Sample from the dynamics for a given input.

% FGPLVM

[mu, var] = gpPosteriorMeanVar(model, X);
y = gsamp(mu, diag(var), 1);
if isfield(model, 'diff') & model.diff
  y = X + y;
end