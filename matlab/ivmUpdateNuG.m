function model = ivmUpdateNuG(model, index)

% IVMUPDATENUG Update nu and g parameters associated with noise model.

% IVM

if nargin < 2
  index = 1:size(model.y, 1);
end

[model.g(index, :), model.nu(index, :)] = ...
    noiseUpdateNuG(model.noise, ...
                   model.mu(index, :), model.varSigma(index, :), ...
                   model.y(index, :));
%/~
if any(model.nu(index, :)< 0) & model.noise.logconcave
  warning('nu less than zero')
end
%~/
