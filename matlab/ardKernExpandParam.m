function kern = ardKernExpandParam(kern, params)

% ARDKERNEXPANDPARAM Create kernel structure from ARD parameters.

% IVM
if kern.linearBound
  params(1:5) = linearBound(params(1:5));
else
  params(1:5) = expBound(params(1:5));
end
kern.inverseWidth = params(1);
kern.rbfVariance = params(2);
kern.whiteVariance = params(3);
kern.biasVariance = params(4);
kern.linearVariance = params(5);

kern.inputScales = sigmoidBound(params(6:end));
