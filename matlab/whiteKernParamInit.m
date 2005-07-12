function kern = whiteKernParamInit(kern)

% WHITEKERNPARAMINIT white noise kernel parameter initialisation.

% KERN

kern.variance = exp(-2);
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = 'negLogLogit';

