function [params, transform] = sqexpKernExtractParam(kern)

% SQEXPKERNEXTRACTPARAM Extract parameters from squared exponential kernel structure.

% IVM

params = [kern.inverseWidth kern.rbfVariance kern.biasVariance ...
          kern.whiteVariance];

