function kern = ardKernParamInit(kern)

% ARDKERNPARAMINIT RBF kerne parameter initialisation.

% IVM

kern.inverseWidth = 1;
kern.rbfVariance = 1;
kern.whiteVariance = 1; 
kern.biasVariance = 1;
kern.linearVariance = 1;
kern.inputScales = ones(1, kern.inputDimension);
kern.nParams = 5 + kern.inputDimension;