% DEMREGRESSION2 Try the IVM for regression.


noiseModel = 'gaussian';
selectionCriterion = 'entropy';
kernelType = 'ARD';
prior = 0;
display = 1;
dVal = 100;
% Sample a regression data-set.
X = randn(1000, 1);
y = [sin(X) + randn(1000, 1)*0.01 -sin(X) + randn(1000, 1)*0.01];

model = ivmRun(X, y, kernelType, noiseModel, selectionCriterion, dVal, prior, display, 100, 4)
model = ivmOptimiseIVM(model, display);

