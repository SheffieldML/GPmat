function kern = linKernExpandParam(kern, params)

% LINKERNEXPANDPARAM Create kernel structure from linear kernel parameters.

% IVM
if kern.linearBound
  kern.variance = linearBound(params(1));
else
  kern.variance = expBound(params(1));
end