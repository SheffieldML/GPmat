load ../../kern/cpp/X.mat
load ../../kern/cpp/y.mat

noiseModel = 'gaussian';
selectionCriterion = 'entropy';

% Just use the rbf ard kernel.
kernelType = {'rbf', 'lin', 'bias', 'white'};

options = ivmOptions;
options.display = 1;
dVal = 10;

model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal);

model = ivmOptimiseIVM(model, options.display);
