function kern = fileKernParamInit(kern)

% FILEKERNPARAMINIT File stored kernel parameter initialisation.

% KERN

kern.variance = 1;
kern.nParams = 1;

kern.transforms.index = [1];
kern.transforms.type = 'negLogLogit';
