% DEMBRENDAN1 Use the GP-LVM to model the Frey face data with back constraints.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'brendan';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('fitc');
options.optimiser = 'conjgrad';
options.back = 'kbr';
options.backOptions = kbrOptions(Y);
options.backOptions.kern = kernCreate(Y, 'rbf');
options.backOptions.kern.inverseWidth = 0.00001;

latentDim = 2;
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


% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'image', [20 28], 1, 0, 1)

% compute the nearest neighbours errors in latent space.
errors = fgplvmNearestNeighbour(model, lbls);
