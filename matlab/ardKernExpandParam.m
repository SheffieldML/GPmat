function kern = ardKernExpandParam(kern, params)

% ARDKERNEXPANDPARAM Create kernel structure from ARD parameters.

% IVM
kern.inverseWidth = params(1);
kern.rbfVariance = params(2);
kern.biasVariance = params(3);
kern.whiteVariance = params(4);
kern.linearVariance = params(5);

kern.inputScales = params(6:end);
