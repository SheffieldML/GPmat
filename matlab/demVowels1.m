% DEMVOWELS1 Model the vowels data with a 2-D FGPLVM using RBF kernel.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'vowels';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
numActive = 200;
latentDim = 2;
model = fgplvmCreate(Y, latentDim, 'fitc', numActive, {'rbf', 'bias', ...
                    'white'}, 'gaussian', 'mlp', 15);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');



% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'vector')

errors = fgplvmNearestNeighbour(model, lbls);
