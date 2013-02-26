function model = gpReversibleDynamicsCreate(q, d, latentVals, options)

% GPREVERSIBLEDYNAMICSCREATE Create a reversible dynamics model.
%
%	Description:
%
%	MODEL = GPREVERSIBLEDYNAMICSCREATE(Q, Q, X, OPTIONS) creates a
%	Gaussian process model for dealing with reversible dynamics in the
%	latent space of a GP-LVM.
%	 Returns:
%	  MODEL - model structure containing the Gaussian process.
%	 Arguments:
%	  Q - the latent space dimension.
%	  Q - the latent space dimension.
%	  X - the latent variables.
%	  OPTIONS - options structure as defined by gpOptions.m.
%	gpReversibleDynamicsLatentGradients, gpReversibleDynamicsSetLatentValues,
%	gpReversibleDynamicsLogLikelihood
%	
%
%	See also
%	GPCREATE, GPREVERSIBLEDYNAMICSOPTIONS, 


%	Copyright (c) 2005 Neil D. Lawrence


if nargin < 4
  options = gpReversibleDynamicsOptions('ftc');
end

diffs = latentVals(2:end, :) - latentVals(1:end-1, :);
X = [latentVals(2:end-1, :) diffs(1:end-1, :)];
y = diffs(2:end, :);
model = gpCreate(2*q, d, X, y, options);
model.learn = 0;
model.type = 'gpReversibleDynamics';

