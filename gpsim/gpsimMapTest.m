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

% SHEFFIELDML

colordef white
load demBarenco1;
B = model.comp{1}.B;
D = model.comp{1}.D;
S = model.comp{1}.S;
clear model

[y, yvar, gene, times, scale, rawExp, rawVar] = gpsimLoadBarencoData;

numGenes = size(gene, 1);
options = gpsimMapOptions(numGenes);
options.intPoints = 17;
options.B = B;
options.D = D;
options.S = S;

nonLinearities = {'linear', 'quadratic', 'negLogLogit', 'exp', 'sigmoid'};
nonLinearities = {'exp'};
for i = 1:length(nonLinearities)
  options.nonLinearity = nonLinearities{i};
  fprintf('Nonlinearity: %s\n', options.nonLinearity);
  model = gpsimMapCreate(numGenes, 1, times, y{1}, yvar{1}, options);
  params = gpsimMapExtractParam(model);
  
  model = gpsimMapExpandParam(model, params);
  
  % Give reasonable values to f.
  f = gsamp(zeros(1, size(model.f, 1)), model.K, 1);
  model = gpsimMapFunctionalExpandParam(model, f);
  
  gradientCheck(params, 'modelObjective', 'modelGradient', model);
  
  model.type = 'gpsimMapFunctional';
  f = gpsimMapFunctionalExtractParam(model);
  gradientCheck(f, 'modelObjective', 'modelGradient', model);
  hessianCheck(f, 'modelObjective', 'modelHessian', model);
end
