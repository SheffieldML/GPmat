function params = mlpardKernExtractParam(kern)

% MLPARDKERNEXTRACTPARAM Extract parameters from multi-layer perceptron ARD kernel structure.

% KERN


params = [kern.weightVariance kern.biasVariance kern.variance kern.inputScales];
