function [params, transform] = sqexpKernExtractParam(kern)

% SQEXPKERNEXTRACTPARAM Extract parameters from squared exponential kernel structure.

% KERN


params = [kern.inverseWidth kern.rbfVariance kern.biasVariance ...
          kern.whiteVariance];

