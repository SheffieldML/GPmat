% DEMBRENDANFGPLVM5 Use the GP-LVM to model the Frey face data with DTCVAR and five latent dimensions..

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'brendan';
experimentNo = 5;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('dtcvar');

latentDim = 5;
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
lvmResultsDynamic(model.type, dataSetName, experimentNo, 'image', [20 28], 1, 0, 1)
