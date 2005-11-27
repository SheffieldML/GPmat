% DEMCLASSIFICATION2 IVM for classification on a data-set sampled from a GP.

% IVM

% Sample a classification data-set.
randn('seed', 1e5)
rand('seed', 1e5)
[X, y ] = mappingLoadData('classificationTwo');
noiseModel = 'probit';
selectionCriterion = 'entropy';
kernelType = {'rbf', 'white'};

options = ivmOptions;
options.display = 2;
dVal = 200;

% Initialise the IVM.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
for i = 1:options.extIters
  % Select the active set.
  model = ivmOptimiseIVM(model, options.display);
  if options.display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Optimise the kernel parameters.
  model = ivmOptimiseKernel(model, options.display, options.kernIters);
end

model = ivmOptimiseIVM(model, options.display);
if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
% Display active points.
model = ivmOptimiseIVM(model, options.display);

% Display the model parameters.
ivmDisplay(model);
