function kern = mlpKernExpandParam(kern, params)

% MLPKERNEXPANDPARAM Create kernel structure from multi-layer perceptron's parameters.

% KERN

%/~
if any(isinf(params))
  warning('params are infinite')
end
%~/

kern.weightVariance = params(1);
kern.biasVariance = params(2);
kern.variance = params(3);
