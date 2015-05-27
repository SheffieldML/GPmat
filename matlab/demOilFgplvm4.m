% DEMOILFGPLVM4 Oil data with deterministic training conditional, and MLP back constraints.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil';
experimentNo = 4;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('dtc');
options.optimiser = 'scg';
options.back = 'mlp';
options.backOptions = mlpOptions;
latentDim = 2;
d = size(Y, 2);

model = fgplvmCreate(latentDim, d, Y, options);

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
lvmScatterPlot(model, lbls);

% compute the nearest neighbours errors in latent space.
errors = lvmNearestNeighbour(model, lbls);
