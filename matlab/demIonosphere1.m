% DEMEP1 Demonstrate Expectation propagation.


randn('seed', 1e6)
rand('seed', 1e6)

% Generate a toy data-set
[X, y] = mapLoadData('ionosphere');
noiseModel = 'probit';
selectionCriterion = 'random';

kernelType = 'rbf';
prior = 0;
display = 0;
dVal = 200;


model = ivm(X, y, kernelType, noiseModel, selectionCriterion, dVal)
model.kern.inverseWidth = 1/(exp(1)*exp(1));
model.kern.variance = exp(4)*exp(4);
model.noise.sigma2 = 1;
%for i = 1:4
%  model = ivmSelectPoints(model, display);
%  model = ivmOptimiseKernel(model, display, 100);
%end
model = ivmSelectPoints(model, display);



