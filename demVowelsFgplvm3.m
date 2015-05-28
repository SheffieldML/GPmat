% DEMVOWELSFGPLVM3 Model the vowels data with a 2-D FGPLVM using RBF kernel and back constraints, but without PCA initialisation.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'vowels';
experimentNo = 3;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('fitc');
options.back = 'mlp';
options.backOptions = mlpOptions;
options.optimiseInitBack = 0;
options.numActive = 200;
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

% Load the results and display dynamically.
lvmResultsDynamic(dataSetName, experimentNo, 'vector')

errors = lvmNearestNeighbour(model, lbls);
