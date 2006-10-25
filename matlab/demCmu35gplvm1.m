% DEMCMU35GPLVM1 Learn a GPLVM on CMU 35 data set.

% FGPLVM

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Get the sequence numbers.
[Y, lbls] = lvmLoadData('cmu35WalkJog');
seq = cumsum(sum(lbls)) - [1:31];

dataSetName = 'cmu35gplvm';
experimentNo = 1;



% load data
[Y, lbls, Ytest, lblstest] = lvmLoadData(dataSetName);
bias = mean(Y);
scale = 1./sqrt(var(Y));
%scale = ones(size(scale));
Y = Y - repmat(bias, size(Y, 1), 1);
Ytest = Ytest - repmat(bias, size(Ytest, 1), 1);
Y = Y.*repmat(scale, size(Y, 1), 1);
Ytest = Ytest.*repmat(scale, size(Ytest, 1), 1);

% Set up model
options = fgplvmOptions('fitc');
options.optimiser = 'conjgrad';
options.back = 'mlp';
options.backOptions = mlpOptions(10);
options.numActive = 200;
options.fixInducing = 1;
options.fixIndices = round(linspace(1, size(Y, 1), options.numActive));
latentDim = 4;

d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);

% Add dynamics model.
optionsDyn = gpOptions('fitc');
optionsDyn.numActive = 200;
options.fixInducing = 1;
diff = 1;
learn = 1;
options.fixIndices = round(linspace(1, size(Y, 1)-length(seq), options.numActive));
model = fgplvmAddDynamics(model, 'gp', optionsDyn, diff, learn, seq);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model', 'scale', 'bias');


