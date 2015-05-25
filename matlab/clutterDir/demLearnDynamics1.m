% DEMLEARNDYNAMICS1 Learn the dynamics associated with the dance.

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'learnDynamics';
experimentNo = 1;

% load data
load dynLearn

% Set up model
numActive = 300;
latentDim = 2;

% Re-create dynamics kernel.
dynKern = kernCreate(model.X, {'rbf', 'white'});
dynKern.comp{1}.inverseWidth = 1;
% This gives signal to noise of 0.1:1e-3 or 100:1.
dynKern.comp{1}.variance = 0.1^2;
dynKern.comp{2}.variance = 1e-3^2;

% use the partially independent training conditional.
model = gpCreate(size(X, 2), size(Y, 2), X, Y, dynKern, 'pitc', numActive);

% Optimise the model.
iters = 1000;
display = 1;

model = gpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

dynKern = model.kern;
save dynKern dynKern
