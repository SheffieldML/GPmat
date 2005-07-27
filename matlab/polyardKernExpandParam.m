function kern = polyardKernExpandParam(kern, params)

% POLYARDKERNEXPANDPARAM Create kernel structure from polynomial ARD's parameters.

% KERN

kern.weightVariance = params(1);
kern.biasVariance = params(2);
kern.variance = params(3);
kern.inputScales = params(4:end);
