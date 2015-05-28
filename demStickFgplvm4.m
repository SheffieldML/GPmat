% DEMSTICKFGPLVM4 Model the stick man using an RBF kernel and 3-D latent space.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'stick';
experimentNo = 4;

% load data
[Y, lbls] = lvmLoadData(dataSetName);

% Set up model
% Train using the full training conditional (i.e. no approximation.)
options = fgplvmOptions('ftc');
latentDim = 3;

d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);
scatter3(model.X(:, 1), model.X(:, 2), model.X(:, 3), 'rx');
% Save the results.
modelWriteResult(model, dataSetName, experimentNo);

