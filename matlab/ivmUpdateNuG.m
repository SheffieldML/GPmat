function model = ivmUpdateNuG(model, index)

% IVMUPDATENUG Update nu and g parameters associated with noise model.

% IVM
% /~[model.nu(index, :), model.g(index, :)] = ...
%     feval([model.noise.type 'NoiseUpdateParams'], ...
%           model.noise, model.mu, ...
%           model.varSigma, model.y, ...
% ~/          index);

if nargin < 2
  index = 1:size(model.y, 1);
end
[model.g(index, :), dlnZ_dvs] = feval([model.noise.type 'NoiseGradVals'], ...
                                      model.noise, ...
                                      model.mu(index, :), ...
                                      model.varSigma(index, :), ...
                                      model.y(index, :));

nu = model.g(index, :).*model.g(index, :) - 2*dlnZ_dvs;
nu(find(abs(nu) < eps)) = eps;
model.nu(index, :) = nu;
if any(model.nu(index, :)< 0)
  warning('nu less than zero')
end

