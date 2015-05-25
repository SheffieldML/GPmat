function classError = ivmGunnarTest(dataSet, dataNum, kernelType, ...
                                    kernParams, noiseParams, beta)

% PPAGUNNARTEST Script for running tests on Gunnar data.

% PPA

% Load the data
HOME = getenv('HOME');
fprintf('Dataset: %s, number %d\n', dataSet, dataNum)
fs = filesep;
baseDir = [HOME filesep 'datasets' filesep 'gunnar' filesep dataSet filesep];
X=load([baseDir dataSet '_train_data_' num2str(dataNum) '.asc']);
y=load([baseDir dataSet '_train_labels_' num2str(dataNum) '.asc']);

% Define the noise model to be used
noiseModel='probit';

% Define the kernel to be used
options = ppaOptions;
options.maxOuterIter = 50;
options.scalarB = 1; % make beta a single scalar.
options.kernIters=0;
options.display = 0;
model=ppa(X, y, noiseModel, kernelType);
model.kern = kernExpandParam(model.kern, kernParams);
model.noise = noiseExpandParam(model.noise, noiseParams);
model.B = repmat(beta, size(model.y));
fprintf('Initial model:\n');
ppaDisplay(model);
model=ppaOptimisePPA(model, options);

fprintf('Final model:\n');
ppaDisplay(model);
fprintf('B: %2.4f\n', model.B(1, 1));
Xtest=load([baseDir dataSet '_test_data_' num2str(dataNum) '.asc']);
ytest=load([baseDir dataSet '_test_labels_' num2str(dataNum) '.asc']);

yPred = ppaOut(model, Xtest);
classError = 1- sum(ytest ==yPred)/length(ytest);

ll = ppaCalculateLogLike2(model);
fprintf('Test Error %2.4f\n', classError);
fprintf('Model likelihood %2.4f\n', ll);

beta = model.B(1, 1);
numIters = model.numIters;
save(['ppa' dataSet num2str(dataNum) 'Test'], 'classError', 'll', ...
     'beta', 'numIters');
fprintf('Data saved ... job finishing ...\n')
     