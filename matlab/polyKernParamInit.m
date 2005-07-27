function kern = polyKernParamInit(kern)

% POLYKERNPARAMINIT Polynomial kernel parameter initialisation.

% KERN

kern.variance = 1;
kern.weightVariance = 1;
kern.biasVariance = 1;
kern.degree = 2;
kern.nParams = 3;

kern.transforms.index = [1 2 3];
kern.transforms.type = 'negLogLogit';
