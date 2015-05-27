% DEMSPGP1DGP4 Do a simple 1-D regression after Snelson & Ghahramani's example.

% GP

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'spgp1d';
experimentNo = 4;

% load data
[X, y] = mapLoadData(dataSetName);

% Set up the model
options = gpOptions('ftc');
options.optimiser = 'conjgrad';

% Make use of correct kernel.
options.kern = kernCreate(X, {'rbf', 'white'});
options.kern.comp{1}.inverseWidth = 20;
options.kern.comp{2}.variance = 0.01;

% Use the full Gaussian process model.
q = size(X, 2);
d = size(y, 2);
model = gpCreate(q, d, X, y, options);

% This updates the kernels.
params = gpExtractParam(model);
model = gpExpandParam(model, params);

% Save results
fileName = modelWriteResult(model, dataSetName, experimentNo);

demSpgp1dPlot
