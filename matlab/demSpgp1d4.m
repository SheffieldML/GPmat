% DEMSPGP1D4 Do a simple 1-D regression after Snelson & Ghahramani's example.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'spgp1d';
experimentNo = 4;

% load data
[X, y] = mapLoadData(dataSetName);

% Set up model
numActive = 9;
latentDim = 2;

% Use the full Gaussian process model.
model = gpCreate(y, X, {'rbf', 'white'}, 'ftc', numActive);
model.kern.comp{1}.inverseWidth = 20;
model.kern.comp{2}.variance = 0.01;

% This updates the kernels.
params = gpExtractParam(model);
model = gpExpandParam(model, params);

% Save results
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


xTest = linspace(-1.5, 1.5, 200)';
[mu, varSigma] = gpPosteriorMeanVar(model, xTest);

figure
plot(X, y, 'r.');
hold on
a = plot(xTest, mu, 'b-');
a = [a plot(xTest, mu+2*sqrt(varSigma), 'b--')];
a = [a plot(xTest, mu-2*sqrt(varSigma), 'b--')];
set(gca, 'ylim', [-1 2])
set(gca, 'xlim', [-1.5 1.5])
set(a, 'linewidth', 2);
zeroAxes
