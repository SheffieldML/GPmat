% DEMBRENDAN1 Model the Frey face data using an RBF kernel and dynamics.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'brendan';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
numActive = 100;
latentDim = 2;

% Train using the full training conditional (i.e. no approximation.)
d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, 'pitc', numActive, {'rbf', 'bias', 'white'}, 'gaussian');

% Add dynamics model.
dynKern = kernCreate(model.X, {'rbf', 'white'});
dynKern.comp{1}.inverseWidth = 1;
% This gives signal to noise of 0.1:1e-3 or 100:1.
dynKern.comp{1}.variance = 0.1^2;
dynKern.comp{2}.variance = 1e-3^2;
model = fgplvmAddDynamics(model, 'gp', dynKern);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'image', [20 28], 1, 0, 1)

