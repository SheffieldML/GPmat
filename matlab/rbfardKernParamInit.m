function kern = rbfardKernParamInit(kern)

% RBFARDKERNPARAMINIT radial basis function ARD kernel parameter initialisation.

% IVM

% This parameter is restricted positive.
kern.inverseWidth = 1;
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 2 + kern.inputDimension;

% Set to 1 to use log(1+exp(a)) to transform parameters in stead of exp(a)
kern.linearBound = 1;