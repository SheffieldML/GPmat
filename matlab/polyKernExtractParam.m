function params = polyKernExtractParam(kern)

% POLYKERNEXTRACTPARAM Extract parameters from polynomial kernel structure.

% KERN

params = [kern.weightVariance kern.biasVariance kern.variance];

%/~
if any(isinf(params))
  warning('params are infinite')
end
%~/
