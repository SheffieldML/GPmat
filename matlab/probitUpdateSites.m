function model = probitUpdateSites(model, index)

% PROBITUPDATESITES Update site parameters for probit model.

% IVM

if nargin < 2
  error('No site index specified');
end

model = probitUpdateParams(model, index);
model = updateSites(model, index);