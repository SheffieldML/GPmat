function kern = rbfKernExpandParam(kern, params)

% RBFKERNEXPANDPARAM Create kernel structure from rbf parameters.

% IVM

if kern.linearBound
  params = linearBound(params);
else
  params = expBound(params);
end
kern.inverseWidth = params(1);
kern.variance = params(2);
