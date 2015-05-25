function model = cmpndUpdateSites(model, index)

% CMPNDUPDATESITES Update site parameters for compound noise model.

% IVM

if nargin < 2
  error('No site index specified');
end

model = cmpndUpdateParams(model, index);
model = updateSites(model, index);