% DEMORDERED2 Run a demonstration of the ordered categories noise model.

% IVM


randn('seed', 1e6)
rand('seed', 1e6)

[X, y] = ivmLoadData('orderedTwo');

noiseModel = 'ordered';
selectionCriterion = 'entropy';
% Use a kernel consisting of an RBF ard kernel, a linear ard kernel and a
% bias term.
kernelType = {'rbfard', 'linard', 'bias', 'white'};
prior = 0;
display = 2;
dVal = 100;

% Initialise the IVM.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

% Constrain the ARD parameters in the RBF and linear kernels to be the same.
model.kern = cmpndTieParameters(model.kern, {[3, 6], [4, 7]});
% Do some plotting
if display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
for i = 1:4
  % Select active set.
  model = ivmOptimiseIVM(model, display);
  if display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Optimise the noise model.
  model = ivmOptimiseNoise(model, prior, display, 100);
  
  % Select active set.
  model = ivmOptimiseIVM(model, display);
  if display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Optimise kernel parameters
  model = ivmOptimiseKernel(model, prior, display, 100);
end
% Select active set.
model = ivmOptimiseIVM(model, display);
if display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
% Display active points
model = ivmOptimiseIVM(model, display);
% Display parameters of end model.
ivmDisplay(model);





