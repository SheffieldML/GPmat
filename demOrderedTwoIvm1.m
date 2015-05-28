% DEMORDEREDTWOIVM1 Run a demonstration of the ordered categories noise model (circular data).

% IVM

randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'orderedTwo';
experimentNo = 1;

% load data
[X, y] = mapLoadData(dataSetName);


% Set up model
options = ivmOptions;
options.noise = 'ordered';
% Learn the noise model for the ordered categorical case.
options.noiseIters = 100;
options.display = 2;
% Use a kernel consisting of an RBF ard kernel, a linear ard kernel and a
% bias term.
options.kern = {'rbfard', 'linard', 'bias', 'white'};

model = ivmCreate(size(X, 1), size(y, 2), X, y, options);

% Constrain the ARD parameters in the RBF and linear kernels to be the same.
model.kern = cmpndTieParameters(model.kern, {[3, 6], [4, 7]});
% Do some plotting
if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
for i = 1:options.extIters
  % Select active set.
  model = ivmOptimiseIvm(model, options.display);
  if options.display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Optimise the noise model.
  model = ivmOptimiseNoise(model, options.display, options.noiseIters);
  
  % Select active set.
  model = ivmOptimiseIvm(model, options.display);
  if options.display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  % Optimise kernel parameters
  model = ivmOptimiseKernel(model, options.display, options.kernIters);
end
% Select active set.
model = ivmOptimiseIvm(model, options.display);
if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
% Display active points
model = ivmOptimiseIvm(model, options.display);
% Display parameters of end model.
ivmDisplay(model);


% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
[kern, noise, ivmInfo] = ivmDeconstruct(model);
save(['dem' capName num2str(experimentNo) '.mat'], ...
     'kern', ...
     'noise', ...
     'ivmInfo');

if exist('printDiagram') & printDiagram
  ivmPrintPlot(model, 'ivmContour', [], [], [], capName, experimentNo);
end
