function [ax, data] = fgplvmDynamicsSample(model, points);

% FGPLVMDYNAMICSSAMPLE Sample a field from the GP.

% FGPLVM

if nargin < 2
  points = 20;
end

% Dynamics samples are only in 2-D at the moment.
if(size(model.X, 2)~=2)
  error(['Latent space should be two-dimensional to sample ' ...
         'dynamics'])
end

fgplvmKernDynamicsSample(model.dynamics.kern, points, model.dynamics.diff);
