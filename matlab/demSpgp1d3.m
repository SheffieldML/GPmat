% DEMSPGP1D3 Do a simple 1-D regression after Snelson & Ghahramani's example.

% GP

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'spgp1d';
experimentNo = 3;

% load data
[X, y] = mapLoadData(dataSetName);

% Set up model
options = gpOptions('pitc');
options.numActive = 9;

% use the deterministic training conditional.
q = size(X, 2);
d = size(y, 2);

model = gpCreate(q, d, X, y, options);
model.X_u = randn(9, 1)*0.25 - 0.75;
params = gpExtractParam(model);
model = gpExpandParam(model, params);

% Optimise the model.
iters = 1000;
display = 1;

model = gpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

demSpgp1dPlot