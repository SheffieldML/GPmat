% DEMTWOCLUSTERS1

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'twoclusters';
experimentNo = 1;

% load data
[X, y] = mapLoadData(dataSetName);

% Set up model
options = gpOptions('fitc');
options.numActive = 100;

% use the deterministic training conditional.
q = size(X, 2);
d = size(y, 2);

model = gpCreate(q, d, X, y, options);
model.X_u = rand(size(model.X_u))*3 - 1.5;
% Optimise the model.
iters = 1000;
display = 1;

model = gpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

demSpgp1dPlot
