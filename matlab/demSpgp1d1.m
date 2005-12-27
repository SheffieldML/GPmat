% DEMSPGP1D1 Do a simple 1-D regression after Snelson & Ghahramani's example.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'spgp1d';
experimentNo = 1;

% load data
[X, y] = mapLoadData(dataSetName);

% Set up model
numActive = 9;
latentDim = 2;

% use the fully independent training conditional.
model = gpCreate(y, X, {'rbf', 'bias', 'white'}, 'dtc', numActive);

model.X_u = randn(9, 1)*0.25 - 0.75;
% Optimise the model.
iters = 1000;
display = 1;

model = gpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


xTest = linspace(-1.5, 1.5, 200)';
[mu, varSigma] = gpPosteriorMeanVar(model, xTest);

figure
plot(X, y, 'r.');
hold on
plot(xTest, mu, 'b-');
plot(xTest, mu+2*sqrt(varSigma), 'b--')
plot(xTest, mu-2*sqrt(varSigma), 'b--')
set(gca, 'ylim', [-1 2])
set(gca, 'xlim', [-1.5 1.5])
zeroAxes