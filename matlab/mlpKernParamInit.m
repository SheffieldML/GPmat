function kern = mlpKernParamInit(kern)

% MLPKERNPARAMINIT multi-layer perceptron kernel parameter initialisation.

% KERN

% KERN


kern.weightVariance = 10;
kern.biasVariance = 10;
kern.variance = 1;
kern.nParams = 3;

kern.transforms.index = [1 2 3];
kern.transforms.type = 'negLogLogit';
% Set to 1 to use log(1+exp(a)) to transform parameters in stead of exp(a)
kern.linearBound = 1;