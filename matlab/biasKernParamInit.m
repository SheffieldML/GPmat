function kern = biasKernParamInit(kern)

% BIASKERNPARAMINIT bias kernel parameter initialisation.

% IVM

kern.variance = 1;
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = 'negLogLogit';
