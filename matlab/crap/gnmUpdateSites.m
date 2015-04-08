function model = gnmUpdateSites(model, index)

% GNMUPDATESITES Update site parameters for gap noise model.

% IVM

if nargin < 2
  error('No site index specified');
end

model = gnmUpdateParams(model, index);
model = updateSites(model, index);