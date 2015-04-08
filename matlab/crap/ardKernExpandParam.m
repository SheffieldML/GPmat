function kern = ardKernExpandParam(params, kern)

% ARDKERNEXPANDPARAM Create kernel structure from ARD parameters.

% IVM

kern.inverseWidth = exp(params(1));
kern.rbfVariance = exp(params(2));
kern.whiteVariance = exp(params(3));
kern.biasVariance = exp(params(4));
kern.linearVariance = exp(params(5));
kern.inputScales = exp(params(6:end));
