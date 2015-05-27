% DEMBRENDANVARGPLVM1 Run variational GPLVM on Brendan face data.

% VARGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'brendan';
experimentNo = 1;
printDiagram = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'white'};
options.numActive = 100; 
options.scale2var1 = 1; % scale data to have variance 1
%options.tieParam = 'tied';  

options.optimiser = 'scg';
latentDim = 30;
d = size(Y, 2);

% create the model
model = vargplvmCreate(latentDim, d, Y, options);
%
model = vargplvmParamInit(model, model.m, model.X); 

% Optimise the model.
iters = 1500;
display = 1;

model = vargplvmOptimise(model, display, iters);

% Save the results.
modelWriteResult(model, dataSetName, experimentNo);

if exist('printDiagram') & printDiagram
  lvmPrintPlot(model, lbls, dataSetName, experimentNo);
end

% Load the results and display dynamically.
lvmResultsDynamic(model.type, dataSetName, experimentNo, 'image', [20 28], 1, 0, 1)
