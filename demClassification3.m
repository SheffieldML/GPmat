% DEMCLASSIFICATION3 IVM for classification on a data-set sampled from a GP with null category.

% IVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'classificationThree';
experimentNo = 3;

% load data
[X, y] = mapLoadData(dataSetName);


% Set up model
options = ivmOptions;
options.display = 2;
options.kern = {'rbf', 'white'};

model = ivmCreate(size(X, 1), size(y, 2), X, y, options);

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







