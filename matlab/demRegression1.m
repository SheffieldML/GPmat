% DEMREGRESSION1 The data-set is sampled from a GP with known parameters.

% IVM

% Sample a regression data-set.
[X, y] = ivmLoadData('regressionOne');


noiseModel = 'gaussian';
selectionCriterion = 'entropy';

% Just use the rbf ard kernel.
kernelType = {'rbfard', 'linard', 'white'};
prior = 0;
display = 2;
dVal = 50;


% Initialise the IVM model.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);
model.kern = cmpndTieParameters(model.kern, {[3, 6], [4, 7]});
if display > 1
  ivm3dPlot(model, 'surf', i);
  colormap bone
  shading interp  
  drawnow
end

for i = 1:4
  % Plot the data.
  % Select the active set.
  model = ivmOptimiseIVM(model, display);
  if display > 1
    ivm3dPlot(model, 'surf', i);
    colormap bone
    shading interp  
    drawnow
  end
  % Optimise kernel parameters.
  model = ivmOptimiseKernel(model, prior, display, 100);

end
model = ivmOptimiseIVM(model, display);
if display > 1
  ivm3dPlot(model, 'surf', i);
  colormap bone
  shading interp  
end
% Show the active points.
model = ivmOptimiseIVM(model, display);

ivmDisplay(model);