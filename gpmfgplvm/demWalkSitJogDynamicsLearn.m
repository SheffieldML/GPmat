% DEMWALKSITJOGDYNAMICSLEARN Learn the stick man dynamics.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'walkSitJog';

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
options = fgplvmOptions('pitc');
options.initX = X; %'isomap';
latentDim = 2;

% Train using the partially independent training conditional.
d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);

% Add dynamics model.
optionsDyn = gpReversibleDynamicsOptions('pitc');
model = fgplvmAddDynamics(model, 'gpReversible', optionsDyn);

model.dynamics.type = 'gp';
display = 1;
iters = 1000;
model.dynamics = gpOptimise(model.dynamics, display, iters);
modle.dynamics.type = 'gpReversibleDynamics';

save(['dem' capName 'DynamicsLearn' '.mat'], 'model');

