function kern = linKernParamInit(kern)

% LINKERNPARAMINIT Linear kernel parameter initialisation.

% KERN

kern.variance = 1;
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = 'negLogLogit';
