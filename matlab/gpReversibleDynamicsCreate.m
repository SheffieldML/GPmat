function model = gpReversibleDynamicsCreate(q, d, latentVals, options)

% GPREVERSIBLEDYNAMICSCREATE Create a reversible dynamics model. 
% FORMAT
% DESC creates a Gaussian process model for dealing with reversible dynamics
% in the latent space of a GP-LVM.
% ARG q : the latent space dimension.
% ARG q : the latent space dimension.
% ARG X : the latent variables.
% ARG options : options structure as defined by gpOptions.m.
% RETURN model : model structure containing the Gaussian process.
%
% SEEALSO : gpCreate, gpReversibleDynamicsOptions,
% gpReversibleDynamicsLatentGradients, gpReversibleDynamicsSetLatentValues,
% gpReversibleDynamicsLogLikelihood
%
% COPYRIGHT : Neil D. Lawrence, 2005

% FGPLVM

if nargin < 4
  options = gpReversibleDynamicsOptions('ftc');
end

diffs = latentVals(2:end, :) - latentVals(1:end-1, :);
X = [latentVals(2:end-1, :) diffs(1:end-1, :)];
y = diffs(2:end, :);
model = gpCreate(2*q, d, X, y, options);
model.learn = 0;
model.type = 'gpReversibleDynamics';

