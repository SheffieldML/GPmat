% DEMSPGP1D2 Do a simple 1-D regression after Snelson & Ghahramani's example.
%
% 

% Copyright (c) 2006 Neil D. Lawrence
% demSpgp1d2.m version 1.3



% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'spgp1d';
experimentNo = 2;

% load data
[X, y] = mapLoadData(dataSetName);

% Set up model
options = gpOptions('fitc');
options.numActive = 9;

% use the deterministic training conditional.
q = size(X, 2);
d = size(y, 2);

model = gpCreate(q, d, X, y, options);

% Optimise the model.
iters = 1000;
display = 1;

model = gpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


demSpgp1dPlot