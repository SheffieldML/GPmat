function kern = linardKernExpandParam(kern, params)

% LINARDKERNEXPANDPARAM Create kernel structure from linear ARD's parameters.

% IVM

if kern.linearBound
  kern.variance = linearBound(params(1));
else
  kern.variance = expBound(params(1));
end
kern.inputScales = sigmoidBound(params(2:end));
