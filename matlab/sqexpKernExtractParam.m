function params = sqexpKernExtractParam(kern)

% SQEXPKERNEXTRACTPARAM Extract parameters from squared exponential kernel structure.

% IVM

params = log([kern.inverseWidth kern.rbfVariance kern.biasVariance kern.whiteVariance]);
