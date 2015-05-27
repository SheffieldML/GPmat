% DEMOILFGPLVM9 Oil data with three dimensions and variational sparse approximation.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil';
experimentNo = 9;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('dtcvar');
options.kern = {'rbf', 'bias', 'whitefixed'};
options.optimiser = 'scg';
latentDim = 3;
d = size(Y, 2);

model = fgplvmCreate(latentDim, d, Y, options);
model.kern.comp{3}.variance = 1e-4;
% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
modelWriteResult(model, dataSetName, experimentNo);

if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, lbls, dataSetName, experimentNo);
end

% Load the results and display them.
lvmScatterPlot(model, lbls,[],[1,2]); % default [1,2], change to plot 
                                      % a different pair of lantent dims

% compute the nearest neighbours errors in latent space.
errors = fgplvmNearestNeighbour(model, lbls);
