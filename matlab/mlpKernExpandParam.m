function kern = mlpKernExpandParam(params, kern)

% MLPKERNEXPANDPARAM Create kernel structure from multi-layer perceptron's parameters.

% IVM

kern.weightVariance = exp(params(1));
kern.biasVariance = exp(params(2));
kern.variance = exp(params(3));
