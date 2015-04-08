% REGRESSIONDEMO Try the IVM for regression.

noiseModel = 'gaussian';
selectionCriterion = 'entropy';
kernelType = 'ARD';
prior = 0;
display = 1;
dVal = 100;
% Sample a regression data-set.
generateRegressionData;

model = ivmRun(X, y, kernelType, noiseModel, selectionCriterion, dVal, prior, display, 100, 4)
model = ivmOptimiseIVM(model, display);

