% DEMROBOTATRIUM1 Wireless Robot data from University of Washington, with back constraints and dynamics.


% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'robotAtrium';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
numActive = 100;
latentDim = 2;

% Train using the full training conditional (i.e. no approximation.)
d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, 'ftc', numActive, {'rbf', 'bias', ...
                    'white'}, 'gaussian', 'mlp', 15);

% Add dynamics model.
model = fgplvmAddDynamics(model, 'gp', {'rbf', 'white'}, 100);
model.dynamics.kern.comp{1}.inverseWidth = 0.2;

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

