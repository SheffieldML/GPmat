function model = ivmUpdateSites(model, index)

% IVMUPDATESITES Update site parameters.

% IVM

[model.m(index, :), model.beta(index, :)] = ...
    noiseUpdateSites(model.noise, ...
                     model.g(index, :), model.nu(index, :), ...
                     model.mu(index, :), model.varSigma(index, :), ...
                     model.y(index, :));


if any(model.beta<0) 
  if model.noise.logconcave
    error('Beta less than zero for log concave model.')
  else
    indices = find(model.beta < 0);
    model.beta(indices) = 0;
    model.m(indices) = 0;
    fprintf('Beta less than zero .... fixing to zero.\n')
  end
end
