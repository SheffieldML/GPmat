function kern = linardKernParamInit(kern)

% LINARDKERNPARAMINIT linear ARD kernel parameter initialisation.

% IVM

% This parameters is restricted positive.
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 1 + kern.inputDimension;

% Set to 1 to use log(1+exp(a)) to transform parameters in stead of exp(a)
kern.linearBound = 1;