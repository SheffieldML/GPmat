function kern = rbfardKernParamInit(kern)

% RBFARDKERNPARAMINIT radial basis function ARD kernel parameter initialisation.

% IVM

kern.inverseWidth = 1;
kern.variance = 1;
kern.inputScales = ones(1, kern.inputDimension);
kern.nParams = 2 + kern.inputDimension;