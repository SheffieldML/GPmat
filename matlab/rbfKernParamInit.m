function kern = rbfKernParamInit(kern)

% RBFKERNPARAMINIT RBF kernel parameter initialisation.

% KERN

kern.inverseWidth = 1;
kern.variance = 1;
kern.nParams = 2;

kern.transforms.index = [1 2];
kern.transforms.type = 'negLogLogit';
