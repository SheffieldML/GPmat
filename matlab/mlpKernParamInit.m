function kern = mlpKernParamInit(kern)

% MLPKERNPARAMINIT multi-layer perceptron kernel parameter initialisation.

% IVM

kern.weightVariance = 1;
kern.biasVariance = 1;
kern.variance = 1;
kern.nParams = 3;