function kern = rbfardKernExpandParam(params, kern)

% RBFARDKERNEXPANDPARAM Create kernel structure from radial basis function ARD's parameters.

% IVM

kern.inverseWidth = exp(params(1));
kern.variance = exp(params(2));
kern.inputScales = sigmoid(params(3:end));
