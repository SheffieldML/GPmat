function model = ivmUpdateSites(model, index)

% IVMUPDATESITES Update site parameters.

% IVM

model.beta(index, :) = 1./(1./model.nu(index, :)-model.varSigma(index, :));
model.m(index, :) = model.mu(index, :) + model.g(index, :)./model.nu(index, :);
