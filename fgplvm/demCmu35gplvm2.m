% DEMCMU35GPLVM2 Learn a GPLVM on CMU 35 data set.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Get the sequence numbers.
[Y, lbls] = lvmLoadData('cmu35WalkJog');
seq = cumsum(sum(lbls)) - [1:31];

dataSetName = 'cmu35gplvm';
experimentNo = 2;

% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
origBias = mean(Y, 1);
origScale = 1./sqrt(var(Y, 1));
Y = Y - repmat(origBias, size(Y, 1), 1);
Ytest = Ytest - repmat(origBias, size(Ytest, 1), 1);
Y = Y.*repmat(origScale, size(Y, 1), 1);
Ytest = Ytest.*repmat(origScale, size(Ytest, 1), 1);

% Set up model
options = fgplvmOptions('fitc');
options.optimiser = 'conjgrad';
options.back = 'mlp';
options.backOptions = mlpOptions(10);
options.numActive = 100;
options.fixInducing = 1;
options.fixIndices = round(linspace(1, size(Y, 1), options.numActive));
latentDim = 3;

d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);

% Add dynamics model.
optionsDyn = gpOptions('fitc');
optionsDyn.numActive = 100;
optionsDyn.fixInducing = 1;
optionsDyn.kern = kernCreate(model.X, {'rbf', 'white'});
optionsDyn.kern.comp{1}.inverseWidth = 0.2;
% This gives signal to noise of 0.1:5e-3 or 20:1.
optionsDyn.kern.comp{1}.variance = 0.01;
optionsDyn.kern.comp{2}.variance = 0.95;
diff = 1;
learn = 1;
optionsDyn.fixIndices = round(linspace(1, size(Y, 1)-length(seq), options.numActive));
model = fgplvmAddDynamics(model, 'gp', optionsDyn, diff, learn, seq);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model', 'origScale', 'origBias');


