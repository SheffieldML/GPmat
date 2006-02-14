function kern = whitefixedKernParamInit(kern)

% WHITEFIXEDKERNPARAMINIT white noise kernel parameter initialisation.

% KERN

kern.variance = exp(-2);
kern.nParams = 0;

%kern.transforms.index = 0;
%kern.transforms.type = 'negLogLogit';