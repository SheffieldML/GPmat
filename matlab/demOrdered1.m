% DEMORDERED1 Run a demonstration of the ordered categories noise model.

% IVM


[X, y] = mappingLoadData('orderedOne');

noiseModel = 'ordered';
selectionCriterion = 'entropy';
% Use a kernel consisting of an RBF ard kernel, a linear ard kernel and a
% bias term.
kernelType = {'rbfard', 'linard', 'bias', 'white'};

options = ivmOptions;
options.noiseIters = 100;
options.display = 2;
dVal = 100;

% Initialise the IVM.
model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

% Constrain the ARD parameters in the RBF and linear kernels to be the same.
model.kern = cmpndTieParameters(model.kern, {[3, 6], [4, 7]});
% Do some plotting
if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
for i = 1:options.extIters
  % Select active set.
  model = ivmOptimiseIVM(model, options.display);
  if options.display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Optimise the noise model.
  model = ivmOptimiseNoise(model, options.display, options.noiseIters);
  
  % Select active set.
  model = ivmOptimiseIVM(model, options.display);
  if options.display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Optimise kernel parameters
  model = ivmOptimiseKernel(model, options.display, options.kernIters);
end
% Select active set.
model = ivmOptimiseIVM(model, options.display);
if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
% Display active points
model = ivmOptimiseIVM(model, options.display);
% Display parameters of end model.
ivmDisplay(model);





