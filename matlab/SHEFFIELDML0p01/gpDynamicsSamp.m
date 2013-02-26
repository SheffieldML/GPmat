function y = gpDynamicsSamp(model, X);

% GPDYNAMICSSAMP Sample from the dynamics for a given input.
%
%	Description:
%
%	GPDYNAMICSSAMP(MODEL, X, X) does a one step ahead sample from the
%	dynamics for a single input location.
%	 Arguments:
%	  MODEL - the dynamics model from which to sample.
%	  X - the input position from which to sample.
%	  X - the new latent position sampled from the dynamics.
%	
%
%	See also
%	GPDYNAMICSCREATE


%	Copyright (c) 2006 Neil D. Lawrence


[mu, var] = gpPosteriorMeanVar(model, X);
y = gsamp(mu, diag(var), 1);
if isfield(model, 'diff') & model.diff
  y = X + y;
end