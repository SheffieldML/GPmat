function params = whiteKernExtractParam(kern)

% WHITEKERNEXTRACTPARAM Extract parameters from white noise kernel structure.

% IVM

params = log(kern.variance);