function kern = sqexpKernParamInit(kern)

% SQEXPKERNPARAMINIT squared exponential kernel parameter initialisation.

% IVM

kern.inverseWidth = 1;
kern.rbfVariance = 1;
kern.whiteVariance = 1; 
kern.biasVariance = 1;
kern.nParams = 4;