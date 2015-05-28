% DEMFOURWALKS1 Model four seperate walsk using an RBF kernel and dynamics.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

dataSetName = 'fourWalks';
experimentNo = 1;

% load data
[Y, lbls] = lvmLoadData(dataSetName);
seq = cumsum(sum(lbls));

% Set up model
% Train using the full training conditional (i.e. no approximation.)
options = fgplvmOptions('ftc');
options.learnScales = 1;
latentDim = 3;

d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);
model.scale = sqrt(var(model.m));
model.m = gpComputeM(model);

% Add dynamics model.
options = gpOptions('ftc');
options.kern = kernCreate(model.X, {'rbf', 'lin', 'white'});
options.kern.comp{1}.inverseWidth = 0.2;
% This gives signal to noise of 0.1:1e-3 or 100:1.
options.kern.comp{1}.variance = 0.009;
options.kern.comp{2}.variance = 0.001;
options.kern.comp{3}.variance = 1e-3^2;
diff = 1;
learn = 0;
model = fgplvmAddDynamics(model, 'gp', options, diff, learn, seq);

% Optimise the model.
iters = 5000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');

demFourWalksReconstruct
