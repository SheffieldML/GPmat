% DEMCLASSIFICATION2 IVM for classification on a data-set sampled from a GP.

% IVM

% Sample a classification data-set.
[X, y ] = ivmLoadData('classificationTwo');
noiseModel = 'probit';
selectionCriterion = 'entropy';
kernelType = {'rbf', 'white'};
prior = 0;
display = 2;
dVal = 200;

% Initialise the IVM.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
if display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
for i = 1:4
  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  if display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Optimise the kernel parameters.
  model = ivmOptimiseKernel(model, prior, display, 100);
end

model = ivmOptimiseIVM(model, display);
if display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
% Display active points.
model = ivmOptimiseIVM(model, display);

% Display the model parameters.
ivmDisplay(model);
