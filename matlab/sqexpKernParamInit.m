function kern = sqexpKernParamInit(kern)

% SQEXPKERNPARAMINIT squared exponential kernel parameter initialisation.

% IVM

kern.inverseWidth = 1;
kern.rbfVariance = 1;
kern.whiteVariance = 1; 
kern.biasVariance = 1;
kern.nParams = 4;

% Set to 1 to use log(1+exp(a)) to transform parameters in stead of exp(a)
kern.linearBound = 1;