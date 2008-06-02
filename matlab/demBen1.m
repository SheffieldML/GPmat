% DEMSPGP1D1 Do a simple 1-D regression after Snelson & Ghahramani's example.

% GP

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'ben';
experimentNo = 1;

load ~/Desktop/FHN_400dps.mat

% load data
X = FHN_400dps_Time';
y = FHN_400dps_NoisyData(1, :)';
% Set up model
options = gpOptions('ftc');
options.kern = {'ratquad', 'white'};

% use the deterministic training conditional.
q = size(X, 2);
d = size(y, 2);

model = gpCreate(q, d, X, y, options);
model.kern.comp{1}.period = 9
%model.X_u = randn(9, 1)*0.25 - 0.75;
%params = gpExtractParam(model);
%model = gpExpandParam(model, params);

% Optimise the model.
iters = 1000;
display = 1;

model = gpOptimise(model, display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');


demSpgp1dPlot