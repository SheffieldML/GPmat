function kern = biasKernParamInit(kern)

% BIASKERNPARAMINIT bias kernel parameter initialisation.

% KERN

% KERN


kern.variance = 1;
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = 'negLogLogit';
