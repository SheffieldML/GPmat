function kern = linardKernParamInit(kern)

% LINARDKERNPARAMINIT linear ARD kernel parameter initialisation.

% IVM

% This parameters is restricted positive.
kern.variance = 2;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.5*ones(1, kern.inputDimension);
kern.nParams = 1 + kern.inputDimension;