function kern = rbfKernParamInit(kern)

% RBFKERNPARAMINIT RBF kernel parameter initialisation.

% IVM

kern.inverseWidth = 1;
kern.variance = 1;
kern.nParams = 2;

% Set to 1 to use log(1+exp(a)) to transform parameters in stead of exp(a)
kern.linearBound = 1;