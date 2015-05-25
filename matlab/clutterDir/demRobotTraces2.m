% DEMROBOTTRACES1 Wireless Robot data from University of Washington, with tailored dynamics.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'robotTraces';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('ftc');
options.kern = {'mlp', 'bias', 'white'};
options.back = 'mlp';
options.backOptions = mlpOptions;
options.optimiseInitBack = 0;
latentDim = 2;


% Train using the full training conditional (i.e. no approximation.)
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



% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'vector')

