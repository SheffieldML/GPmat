function kern = rbfKernExpandParam(params, kern)

% RBFKERNEXPANDPARAM Create kernel structure from rbf parameters.

% IVM

kern.inverseWidth = exp(params(1));
kern.variance = exp(params(2));
