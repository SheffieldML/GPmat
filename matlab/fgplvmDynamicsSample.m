function [ax, data] = fgplvmDynamicsSample(model, points);

% FGPLVMDYNAMICSSAMPLE Sample a field from the GP.

% FGPLVM

% Dynamics samples are only in 2-D at the moment.
if(size(model.X, 2)~=2)
  error(['Latent space should be two-dimensional to sample ' ...
         'dynamics'])
end

if nargin < 2
  fgplvmKernDynamicsSample(model.dynamics.kern);
else
  fgplvmKernDynamicsSample(model.dynamics.kern, points);
end