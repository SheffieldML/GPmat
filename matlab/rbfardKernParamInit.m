function kern = rbfardKernParamInit(kern)

% RBFARDKERNPARAMINIT radial basis function ARD kernel parameter initialisation.

% IVM

% This parameter is restricted positive.
kern.inverseWidth = 2;
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.5*ones(1, kern.inputDimension);
kern.nParams = 2 + kern.inputDimension;