function model = ivmUpdateSites(model, index)

% IVMUPDATESITES Update site parameters.

% IVM

[model.m(index, :), model.beta(index, :)] = ...
    noiseUpdateSites(model.noise, ...
                     model.g(index, :), model.nu(index, :), ...
                     model.mu(index, :), model.varSigma(index, :), ...
                     model.y(index, :));


if any(model.beta<0)
  warning('Beta less than zero')
end
