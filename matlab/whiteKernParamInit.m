function kern = whiteKernParamInit(kern)

% WHITEKERNPARAMINIT white noise kernel parameter initialisation.

% IVM

kern.variance = 1;
kern.nParams = 1;

% Set to 1 to use log(1+exp(a)) to transform parameters in stead of exp(a)
kern.linearBound = 1;