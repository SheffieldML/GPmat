function kern = whiteKernExpandParam(kern, params)

% WHITEKERNEXPANDPARAM Create kernel structure from white noise's parameters.

% IVM

kern.variance = params(1);