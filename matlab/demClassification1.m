% DEMCLASSIFICATION1 Test IVM code on a toy feature selection

% IVM

randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
X = [randn(100,2)-[zeros(100, 1) 6*ones(100, 1)]; randn(100,2)+[zeros(100, 1) 6*ones(100, 1)]; randn(100, 2)];
y = [ones(200, 1); -ones(100, 1)];

% The probit is a classification noise model.
noiseModel = 'probit'; 
selectionCriterion = 'entropy';

% Use a combination of an MLP and linear ARD kernel.
kernelType = {'mlpard', 'linard', 'white'};
dVal = 100;

options = ivmOptions;
options.display = 2;

% Initialise the model.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

% Constrain the ARD parameters in the MLP and linear kernels to be the same.
model.kern = cmpndTieParameters(model.kern, {[4, 7], [5, 8]});

if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
for i = 1:options.extIters;

  % Select the active set.
  model = ivmOptimiseIVM(model, options.display);
  % Plot the data.
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
% display active points.
model = ivmOptimiseIVM(model, options.display);

% Display the final model.
ivmDisplay(model);





