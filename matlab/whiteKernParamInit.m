function kern = whiteKernParamInit(kern)

% WHITEKERNPARAMINIT white noise kernel parameter initialisation.

% IVM

kern.variance = 1;
kern.nParams = 1;

kern.transforms.index = 1;
kern.transforms.type = 'negLogLogit';

