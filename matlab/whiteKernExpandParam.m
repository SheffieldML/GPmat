function kern = whiteKernExpandParam(params, kern)

% WHITEKERNEXPANDPARAM Create kernel structure from white noise's parameters.

% IVM

kern.variance = exp(params(1));