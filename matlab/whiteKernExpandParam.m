function kern = whiteKernExpandParam(kern, params)

% WHITEKERNEXPANDPARAM Create kernel structure from white noise's parameters.

% IVM

if kern.linearBound
  kern.variance = linearBound(params(1));
else
  kern.variance = expBound(params(1));
end