% DEMREGRESSION1 The data-set is sampled from a GP with known parameters.

% IVM

noiseModel = 'gaussian';
selectionCriterion = 'entropy';

% Just use the rbf ard kernel.
kernelType = {'rbfard', 'white'};
prior = 0;
display = 2;
dVal = 50;

% Sample a regression data-set.
generateRegressionData;

% Initialise the IVM model.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

for i = 1:4
  % Plot the data.
  if display > 1
    ivm3dPlot(model, 'surf', i);
    colormap bone
    shading interp  
  end
  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  % Optimise kernel parameters.
  model = ivmOptimiseKernel(model, prior, display, 100);

end
if display > 1
  ivm3dPlot(model, 'surf', i);
  colormap bone
  shading interp  
end
model = ivmOptimiseIVM(model, display);
ivmDisplay(model);