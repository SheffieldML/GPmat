function model = orderedUpdateSites(model, index)

% ORDEREDUPDATESITES Update site parameters for ordered model.

% IVM

if nargin < 2
  error('No site index specified');
end

model = orderedUpdateParams(model, index);
model = updateSites(model, index);