function params = mlpKernExtractParam(kern)

% MLPKERNEXTRACTPARAM Extract parameters from multi-layer perceptron kernel structure.

% IVM
if kern.linearBound
  params = log(exp([kern.weightVariance kern.biasVariance kern.variance])-1);
else
  params = log([kern.weightVariance kern.biasVariance kern.variance]);
end