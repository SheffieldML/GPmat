function kern = rbfardKernExpandParam(kern, params)

% RBFARDKERNEXPANDPARAM Create kernel structure from radial basis function ARD's parameters.

% IVM
kern.inverseWidth = params(1);
kern.variance = params(2);
kern.inputScales = params(3:end);
