function kern = rbfKernParamInit(kern)

% RBFKERNPARAMINIT RBF kernel parameter initialisation.

% IVM

kern.inverseWidth = 1;
kern.variance = 1;
kern.nParams = 2;