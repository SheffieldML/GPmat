function kern = ardKernParamInit(kern)

% ARDKERNPARAMINIT RBF kerne parameter initialisation.

% IVM

kern.inverseWidth = 1;
kern.rbfVariance = 1;
kern.whiteVariance = 1; 
kern.biasVariance = 1;
kern.linearVariance = 1;
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 5 + kern.inputDimension;

% Set to 1 to use log(1+exp(a)) to transform parameters in stead of exp(a)
kern.linearBound = 1;