function model = heavisideUpdateSites(model, index)

% HEAVISIDEUPDATESITES Update site parameters for heaviside model.

% IVM

if nargin < 2
  error('No site index specified');
end

model = heavisideUpdateParams(model, index);
model = updateSites(model, index);
