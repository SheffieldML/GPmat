function kern = linardKernExpandParam(kern, params)

% LINARDKERNEXPANDPARAM Create kernel structure from linear ARD's parameters.

% IVM

kern.variance = params(1);
kern.inputScales = params(2:end);
