function params = mlpKernExtractParam(kern)

% MLPKERNEXTRACTPARAM Extract parameters from multi-layer perceptron kernel structure.

% KERN


params = [kern.weightVariance kern.biasVariance kern.variance];

%/~
if any(isinf(params))
  warning('params are infinite')
end
%~/
