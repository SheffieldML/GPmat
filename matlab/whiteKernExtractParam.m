function params = whiteKernExtractParam(kern)

% WHITEKERNEXTRACTPARAM Extract parameters from white noise kernel structure.

% IVM

if kern.linearBound
  params = log(exp(kern.variance)-1);
else
  params = log(kern.variance);
end