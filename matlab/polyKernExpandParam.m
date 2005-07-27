function kern = polyKernExpandParam(kern, params)

% POLYKERNEXPANDPARAM Create kernel structure from polynomial kernel parameters.
% KERN

%/~
if any(isinf(params))
  warning('params are infinite')
end
%~/

kern.weightVariance = params(1);
kern.biasVariance = params(2);
kern.variance = params(3);
