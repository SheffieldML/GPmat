% DEMCLASSIFICATION2 IVM for classification on a data-set sampled from a GP.

% IVM

importTool('prior');
importTool('kern');
importTool('noise');
importTool('optimi');

noiseModel = 'probit';
selectionCriterion = 'entropy';
kernelType = {'rbf', 'white'};
prior = 0;
display = 2;
dVal = 200;
% Sample a classification data-set.
generateClassificationData;

% Initialise the IVM.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
for i = 1:4
  if display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  % Optimise the kernel parameters.
  model = ivmOptimiseKernel(model, prior, display, 100);
end

model = ivmOptimiseIVM(model, display);
% Display the model parameters.
ivmDisplay(model);
