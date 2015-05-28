% DEMSTICKFGPLVM3 Model the stick man using an RBF kernel and RBF kernel based back constraints.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'stick';
experimentNo = 3;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
% Train using the full training conditional (i.e. no approximation.)
options = fgplvmOptions('ftc');
latentDim = 2;
options.back = 'kbr';
options.backOptions = kbrOptions(Y);
options.backOptions.kern = kernCreate(Y, 'rbf');
options.backOptions.kern.inverseWidth = 0.0001;
options.optimiseInitBack = true;

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


% load connectivity matrix
[void, connect] = mocapLoadTextData('run1');
% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'stick', connect)

