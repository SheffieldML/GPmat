% DEMEP1 Demonstrate Expectation propagation on a toy data set..

% IVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'classificationOne';
experimentNo = 2;

% load data
[X, y] = mapLoadData(dataSetName);




options = ivmOptions;
options.kern = {'mlp', 'white'};
options.display = 2;
options.numActive = 100;

model = ivmCreate(size(X, 1), size(y, 2), X, y, options);

if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
for i= 1:4
  model = ivmSelectPoints(model, options.display);
  % Plot the data.
  if options.display > 1
    ivm3dPlot(model, 'ivmContour', i);
  end
  model = ivmOptimiseKernel(model, options.display, options.kernIters);
end
if options.display > 1
  ivm3dPlot(model, 'ivmContour', i);
end
% display active points.
model = ivmSelectPoints(model, options.display);



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


