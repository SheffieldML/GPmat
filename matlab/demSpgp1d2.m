% DEMSPGP1D2 Do a simple 1-D regression after Snelson & Ghahramani's example.

% FGPLVM

% Fix seeds
randn('seed', 2e5);
rand('seed', 2e5);
seedVal = 2e5;
dataSetName = 'spgp1d';
experimentNo = 2;

% load data
[X, y] = mapLoadData(dataSetName, seedVal);

% Set up model
options = gpOptions('fitc');
options.kern = 'rbf';
options.numActive = 9;
options.optimiser = 'conjgrad';
% use the deterministic training conditional.
q = size(X, 2);
d = size(y, 2);

model = gpCreate(q, d, X, y, options);
model.betaTransform = 'exp';
model.kern.transforms.type = 'exp';
% Optimise the model.
iters = 1000;
display = 1;
model.beta = 4/var(y);
model.kern.variance = var(y);
model.kern.inverseWidth = 1./((-min(X)+max(X))'/2).^2

model = gpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


demSpgp1dPlot