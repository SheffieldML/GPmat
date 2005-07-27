function params = polyardKernExtractParam(kern)

% POLYARDKERNEXTRACTPARAM Extract parameters from multi-layer perceptron ARD kernel structure.

% KERN

params = [kern.weightVariance kern.biasVariance kern.variance kern.inputScales];
