function kern = sqexpKernParamInit(kern)

% SQEXPKERNPARAMINIT squared exponential kernel parameter initialisation.

% KERN

% KERN


kern.inverseWidth = 1;
kern.rbfVariance = 1;
kern.whiteVariance = 1; 
kern.biasVariance = 1;
kern.nParams = 4;


kern.transforms(1).index = [1 2 3 4];
kern.transforms(1).type = 'negLogLogit';
