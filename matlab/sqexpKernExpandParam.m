function kern = sqexpKernExpandParam(kern, params)

% SQEXPKERNEXPANDPARAM Create kernel structure from squared exponential's parameters.

% KERN

kern.inverseWidth = params(1);
kern.rbfVariance = params(2);
kern.biasVariance = params(3);
kern.whiteVariance = params(4);
