% DEMDANCE2 Model the dance man using an RBF kernel and dynamics.


% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'dance';
experimentNo = 2;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
numActive = 300;
latentDim = 2;

load demDance1.mat

% Load learnt dynamics model, from learning on demDance1 dynamics.
load dynKern
model.type = 'fgplvm';
model = fgplvmAddDynamics(model, 'gp', dynKern);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

% load connectivity matrix
[void, connect] = mocapLoadTextData('dancemotion');
% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'stick', connect)
