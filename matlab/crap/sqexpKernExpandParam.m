function kern = sqexpKernExpandParam(params, kern)

% SQEXPKERNEXPANDPARAM Create kernel structure from squared exponential's parameters.

% IVM

kern.inverseWidth = exp(params(1));
kern.rbfVariance = exp(params(2));
kern.whiteVariance = exp(params(3));
kern.biasVariance = exp(params(4));
