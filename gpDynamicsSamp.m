function y = gpDynamicsSamp(model, X);

% GPDYNAMICSSAMP Sample from the dynamics for a given input.
% FORMAT
% DESC does a one step ahead sample from the dynamics for a single
% input location.
% ARG model : the dynamics model from which to sample.
% ARG X : the input position from which to sample.
% ARG X : the new latent position sampled from the dynamics.
% 
% SEEALSO : gpDynamicsCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

[mu, var] = gpPosteriorMeanVar(model, X);
y = gsamp(mu, diag(var), 1);
if isfield(model, 'diff') & model.diff
  y = X + y;
end
