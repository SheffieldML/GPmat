function model = gaussianUpdateSites(model, index)

% GAUSSIANUPDATESITES Update sites for gaussian noise model.

% IVM

if nargin < 2
  error('No site index specified');
end

model = gaussianUpdateParams(model, index);
model.beta(index, :) = 1./model.noise.sigma2;
model.m(index, :) = model.y(index, :)-model.noise.bias;
% model.beta and model.m are already known and set