function params = mlpKernExtractParam(kern)

% MLPKERNEXTRACTPARAM Extract parameters from multi-layer perceptron kernel structure.

% IVM

params = log([kern.weightVariance kern.biasVariance kern.variance]);