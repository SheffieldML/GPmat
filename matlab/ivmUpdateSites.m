function model = ivmUpdateSites(model, index)

% IVMUPDATESITES Update site parameters.

% IVM
  
model.beta(index, :) = model.nu(index, :) ...
    ./(1 - model.nu(index, :).*model.varSigma(index, :));
model.m(index, :) = model.mu(index, :) + model.g(index, :)./model.nu(index, :);

if any(model.beta<0)
  warning('Beta less than zero')
end
