function kern = sqexpKernExpandParam(kern, params)

% SQEXPKERNEXPANDPARAM Create kernel structure from squared exponential's parameters.

% IVM
if kern.linearBound
  params = linearBound(params);
else
  params = expBound(params);
end
kern.inverseWidth = params(1);
kern.rbfVariance = params(2);
kern.whiteVariance = params(3);
kern.biasVariance = params(4);
