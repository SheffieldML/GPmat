function kern = mlpKernParamInit(kern)

% MLPKERNPARAMINIT multi-layer perceptron kernel parameter initialisation.

% KERN

kern.weightVariance = 10;
kern.biasVariance = 10;
kern.variance = 1;
kern.nParams = 3;

kern.transforms.index = [1 2 3];
kern.transforms.type = 'negLogLogit';
