function model = gpsimMapTest

% GPSIMMAPTEST Tests the gradients of the GPSIMMAP model.
% FORMAT
% DESC tests the gradients and log likelihoods of the MAP
% approximation for the GPSIM model.
% RETURN model : the model that was used for testing.
% 
% SEEALSO : gpsimMapCreate
% 
% COPYRIGHT : Neil D. Lawrence, 2006

% GPSIM

numGenes = 5;
times = [0:2:12]';
y = exp(randn(size(times, 1), numGenes));
yvar = exp(randn(size(times, 1), numGenes));
options = gpsimMapOptions(numGenes);
options.intPoints = 17;
nonLinearities = {'linear', 'exp'};
for i = 1:length(nonLinearities)
  options.nonLinearity = nonLinearities{i};
  fprintf('Nonlinearity: %s\n', options.nonLinearity);
  model = gpsimMapCreate(numGenes, 1, times, y, yvar, options);
  params = gpsimMapExtractParam(model);
  params = 10*rand(size(params));
  model = gpsimMapExpandParam(model, params);
  
  % Give reasonable values to f.
  model = gpsimMapFunctionalExpandParam(model, ...
                                        gsamp(zeros(1, size(model.f, 1)), ...
                                              model.K, 1));
  
  
  gradientCheck(params, 'modelObjective', 'modelGradient', model);
  
  model.type = 'gpsimMapFunctional';
  f = gpsimMapFunctionalExtractParam(model);
  gradientCheck(f, 'modelObjective', 'modelGradient', model);
  hessianCheck(f, 'modelObjective', 'modelHessian', model);
end