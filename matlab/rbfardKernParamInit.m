function kern = rbfardKernParamInit(kern)

% RBFARDKERNPARAMINIT radial basis function ARD kernel parameter initialisation.

% KERN

% This parameter is restricted positive.
kern.inverseWidth = 2/kern.inputDimension;
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 2 + kern.inputDimension;

kern.transforms(1).index = [1 2];
kern.transforms(1).type = 'negLogLogit';
kern.transforms(2).index = [3:kern.nParams];
kern.transforms(2).type = 'sigmoid';
