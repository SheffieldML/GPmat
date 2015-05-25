% DEMOILEZFT Model the oil data with the EZFT algorithm.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil100';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);
% Set up model
numActive = 50;
latentDim = 2;
model = fgplvmCreate(Y, latentDim, numActive, {'rbf', 'bias', 'white'}, 'gaussian');

% Plot the intialisation.
symbols = getSymbols(3);
figure, hold on
for i = 1:size(model.X, 1)
  labelNo = find(lbls(i, :));
  plot(model.X(i, 1), model.X(i, 2), symbols{labelNo})
end

% Optimise the model.
iters = 500;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');



% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'vector')

errors = fgplvmNearestNeighbour(model, lbls);

% Load the results and display statically.
% gplvmEzftResultsStatic(dataSetName, experimentNo, 'vector')

% Load the results and display as scatter plot
% gplvmEzftResultsStatic(dataSetName, experimentNo, 'none')
