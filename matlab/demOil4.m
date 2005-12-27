% DEMOIL5 Oil data with deterministic training conditional, and MLP back constraints.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'oil';
experimentNo = 4;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
numActive = 100;
latentDim = 2;

% use the deterministic training conditional back constrained by mlp.
model = fgplvmCreate(Y, latentDim, 'dtc', numActive, {'rbf', 'bias', ...
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

% compute the nearest neighbours errors in latent space.
errors = fgplvmNearestNeighbour(model, lbls);
