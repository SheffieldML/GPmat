function kern = mlpKernExpandParam(kern, params)

% MLPKERNEXPANDPARAM Create kernel structure from multi-layer perceptron's parameters.

% IVM
if kern.linearBound
  params = linearBound(params);
else
  params = expBound(params);
end

kern.weightVariance = params(1);
kern.biasVariance = params(2);
kern.variance = params(3);
