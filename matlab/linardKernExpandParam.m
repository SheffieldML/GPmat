function kern = linardKernExpandParam(kern, params)

% LINARDKERNEXPANDPARAM Create kernel structure from linear ARD's parameters.

% KERN


kern.variance = params(1);
kern.inputScales = params(2:end);
