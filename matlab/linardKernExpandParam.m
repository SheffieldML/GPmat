function kern = linardKernExpandParam(params, kern)

% LINARDKERNEXPANDPARAM Create kernel structure from linear ARD's parameters.

% IVM


kern.variance = exp(params(1));
kern.inputScales = sigmoid(params(2:end));
