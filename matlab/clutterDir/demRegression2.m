% DEMREGRESSION2 The data-set is sampled from a GP with known parameters.

% IVM

% Sample a regression data-set.
[X, y] = mappingLoadData('regressionTwo');

noiseModel = 'gaussian';
selectionCriterion = 'entropy';

% Just use the rbf ard kernel.
kernelType = {'rbfard', 'linard', 'white'};

options = ivmOptions;
options.display = 2;
dVal = 100;

i = 0;
% Initialise the IVM model.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
model.kern = cmpndTieParameters(model.kern, {[3, 6], [4, 7]});

if options.display > 1
  [h1, h2] = ivm3dPlot(model, 'mesh', i);
  drawnow
end

for i = 1:options.extIters
  % Plot the data.
  % Select the active set.
  model = ivmOptimiseIVM(model, options.display);
  if options.display > 1
    delete(h2)
    [h1, h2] = ivm3dPlot(model, 'mesh', i);
    drawnow
  end
  % Optimise kernel parameters.
  model = ivmOptimiseKernel(model, options.display, options.kernIters);

end
model = ivmOptimiseIVM(model, options.display);
if options.display > 1
  delete(h2)
  [h1, h2] = ivm3dPlot(model, 'mesh', i);
end
% Show the active points.
model = ivmOptimiseIVM(model, options.display);

ivmDisplay(model);