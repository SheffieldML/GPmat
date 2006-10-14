% DEMROBOTTRACES1 Wireless Robot data from University of Washington, with tailored dynamics.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'robotTraces';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('ftc');
options.isSpherical = 0;
options.isMissingData = 1;
latentDim = 2;


% Train using the full training conditional (i.e. no approximation.)
d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);

% Add dynamics model.
optionsDyn = gpReversibleDynamicsOptions('ftc');
options.kern.comp{1}.comp{1}.inverseWidth = 1;
model = fgplvmAddDynamics(model, 'gpReversible', optionsDyn);

% Optimise the model.
iters = 10000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');



% Load the results and display dynamically.
fgplvmResultsDynamic(dataSetName, experimentNo, 'vector')

