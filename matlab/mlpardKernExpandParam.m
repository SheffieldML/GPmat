function kern = mlpardKernExpandParam(kern, params)

% MLPARDKERNEXPANDPARAM Create kernel structure from multi-layer perceptron ARD's parameters.

% KERN

% KERN


kern.weightVariance = params(1);
kern.biasVariance = params(2);
kern.variance = params(3);
kern.inputScales = params(4:end);
