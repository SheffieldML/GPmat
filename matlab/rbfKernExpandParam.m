function kern = rbfKernExpandParam(kern, params)

% RBFKERNEXPANDPARAM Create kernel structure from rbf parameters.

% IVM

kern.inverseWidth = params(1);
kern.variance = params(2);
