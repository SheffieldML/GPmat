function kern = linKernExpandParam(params, kern)

% LINKERNEXPANDPARAM Create kernel structure from linear kernel parameters.

% IVM

kern.variance = exp(params(1));