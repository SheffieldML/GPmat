function kern = mlpardKernParamInit(kern)

% MLPARDKERNPARAMINIT multi-layer perceptron ARD kernel parameter initialisation.

% KERN

kern.weightVariance = 10;
kern.biasVariance = 10;
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 3 + kern.inputDimension;

kern.transforms(1).index = [1 2 3];
kern.transforms(1).type = 'negLogLogit';
kern.transforms(2).index = [4:kern.nParams];
kern.transforms(2).type = 'sigmoid';

