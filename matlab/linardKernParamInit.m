function kern = linardKernParamInit(kern)

% LINARDKERNPARAMINIT linear ARD kernel parameter initialisation.

% IVM

kern.variance = 1;
kern.inputScales = ones(1, kern.inputDimension);
kern.nParams = 1 + kern.inputDimension;