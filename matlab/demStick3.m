% DEMSTICK3 Model the stick man using an RBF kernel and RBF kernel based back constraints.

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

d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

if exist('printDiagram') & printDiagram
  fgplvmPrintPlot(model, lbls, capName, experimentNo);
end

% load connectivity matrix
[void, connect] = mocapLoadTextData('run1');
% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'stick', connect)

